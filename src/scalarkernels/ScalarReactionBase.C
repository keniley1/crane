#include "ScalarReactionBase.h"

//registerMooseObject("CraneApp", ScalarReactionBase);

template <>
InputParameters
validParams<ScalarReactionBase>()
{
  InputParameters params = validParams<ODEKernel>();
  params.addRequiredCoupledVar("reactants",
                               "All of the reactants involved in this reaction. If this is a "
                               "townsend reaction, electrons are excluded.");
  params.addRequiredParam<std::vector<Real>>(
      "coefficients",
      "The number of each reactant that participates in this reaction. For example, in the "
      "ionization reaction, e + Ar -> 2e + Ar+, the coefficient of both e and Ar is 1 because only "
      "one of each appears on the reactant side of the equation.");
  params.addRequiredParam<std::string>("reaction", "Stores the full reaction equation.");
  params.addParam<std::string>(
      "number",
      "",
      "The reaction number. Optional, just for material property naming purposes. If a single "
      "reaction has multiple different rate coefficients (frequently the case when multiple "
      "species are lumped together to simplify a reaction network), this will prevent the same "
      "material property from being declared multiple times.");
  params.addClassDescription("Base class for generic reactions. Multiplies the reactants of the "
                             "specified reaction together to be supplied to the child class.");
  return params;
}

ScalarReactionBase::ScalarReactionBase(const InputParameters & parameters)
  : ODEKernel(parameters),
  /*
    _v_var(isCoupledScalar("v") ? coupledScalar("v") : -1),
    _v(coupledScalarValue("v")),
    _rate_coefficient(coupledScalarValue("rate_coefficient")),
    _stoichiometric_coeff(getParam<Real>("coefficient")),
    _v_eq_u(getParam<bool>("v_eq_u")),
    _rate_constant_equation(getParam<bool>("rate_constant_equation"))
    */
    _coefficients(getParam<std::vector<Real>>("coefficients"))
{
  int dval = -1;
  _num_reactants = coupledScalarComponents("reactants");

  /*
  if (_townsend && _num_reactants > 2)
    mooseError("Only two body Townsend reactions are supported. Please make sure that the reaction "
               "is written correctly, and if so, set townsend = false for this reaction.");
               */

  _reactants.resize(_num_reactants);
  _reactants_var.resize(_num_reactants);

  for (unsigned int i = 0; i < _num_reactants; ++i)
  {
    // If this is a townsend reaction, check to make sure electrons are
    // not included in the reactant list. Those are multiplied in later.
    // If included, it is ignored and that element is erased from
    // the reactant vector.
    // Foolproof! (Hopefully...)
    /*
    if (_townsend && getScalarVar("reactants", i)->name() == getParam<std::string>("electron_name"))
    {
      dval = i;
      continue;
    }
    */

    _reactants[i] = &coupledScalarValue("reactants", i);
    _reactants_var[i] = coupledScalar("reactants", i);
  }

  // Delete the element corresponding to electrons from reactant list if this
  // is labeled as a Townsend reaction.
  if (dval >= 0)
  {
    if (_coefficients.size() == _reactants.size())
      _coefficients.erase(_coefficients.begin() + dval);

    _reactants.erase(_reactants.begin() + dval);
    _num_reactants += -1;
  }

  if (_num_reactants != _coefficients.size())
    mooseError("The number of reactants and coefficients is not equal!");
}

Real
ScalarReactionBase::multiplyReactants()
{
  _val = 1;
  for (unsigned int i = 0; i < _num_reactants; ++i)
  {
    _val *= std::pow((*_reactants[i])[_i], _coefficients[i]);
  }

  return _val;
}

Real
ScalarReactionBase::multiplyReactantsJacobian(unsigned int ivar)
{
  _val = 1;
  bool count = false;
  for (unsigned int i = 0; i < _num_reactants; ++i)
  {
    if (_reactants_var[i] == ivar)
    {
      _val *= std::pow((*_reactants[i])[_i], _coefficients[i] - 1);
      count = true;
    }
    else
      _val *= std::pow((*_reactants[i])[_i], _coefficients[i]);
  }

  // The 'count' variable stores a boolean ensuring that the _u variable
  // is actually included in the reactants.
  // If it is, count == true and the jacobian is computed.
  // If not, count == false and the jacobian is zero. 

  if (count)
    return _val;
  else
    return 0;
}

Real
ScalarReactionBase::multiplyReactantsDerivative(unsigned int ivar)
{
  _val = 1;
  bool count = false;
  for (unsigned int i = 0; i < _num_reactants; ++i)
  {
    if (_reactants_var[i] == ivar)
    {
      _val *= _coefficients[i] * std::pow((*_reactants[i])[_i], _coefficients[i] - 1);
      count = true;
    }
    else
      _val *= std::pow((*_reactants[i])[_i], _coefficients[i]);
  }

  // The 'count' variable stores a boolean ensuring that the active variable
  // is actually included in the reactants.
  // If it is, count == true and the jacobian is computed.
  // If not, count == false and the jacobian is zero. 

  if (count)
    return _val;
  else
    return 0;
}

/*
Real
ScalarReactionBase::computeQpResidual()
{
  return -_stoichiometric_coeff * _rate_coefficient[_i] * _v[_i];
}

Real
ScalarReactionBase::computeQpJacobian()
{
  Real power, eq_u_mult, gas_mult;

  power = 0.0;
  eq_u_mult = 1.0;
  gas_mult = 1.0;

  if (isCoupledScalar("v") && _v_eq_u)
  {
    power += 1.0;
    eq_u_mult = _v[_i];
  }
  else
    gas_mult *= _v[_i];

  return -_stoichiometric_coeff * _rate_coefficient[_i] * gas_mult * power *
         std::pow(eq_u_mult, power - 1);
}

Real
ScalarReactionBase::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real power, deriv_factor, other_factor;
  power = 0;
  other_factor = 1;
  deriv_factor = 1;

  if (jvar == _v_var && !_v_eq_u)
  {
    power += 1;
    deriv_factor = _v[_i];
  }
  else
    other_factor *= _v[_i];

  return -_stoichiometric_coeff * _rate_coefficient[_i] * other_factor * power *
         std::pow(deriv_factor, power - 1);
}
*/
