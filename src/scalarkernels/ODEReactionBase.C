#include "ODEReactionBase.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"

// registerMooseObject("CraneApp", ODEReactionBase);

defineLegacyParams(ODEReactionBase);

InputParameters
ODEReactionBase::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addRequiredCoupledVar("reactants", "All coupled variables in the syste1.");
  params.addRequiredParam<Real>("coefficient", "The stoichiometric coefficient.");
  params.addParam<bool>(
      "rate_constant_equation", false, "True if rate constant is provided by equation.");
  return params;
}

ODEReactionBase::ODEReactionBase(const InputParameters & parameters)
  : ODEKernel(parameters),
    _num_reactants(coupledScalarComponents("reactants")),
    _reactant_var(coupledScalar("reactants")),
    _reactants(_num_reactants),
    _reactant_num(_num_reactants),
    //_reactants(coupledScalarValue("reactants")),
    _stoichiometric_coeff(getParam<Real>("coefficient"))
{
  for (unsigned int i = 0; i < _num_reactants; ++i)
  {
    _reactants[i] = &coupledScalarValue("reactants", i);
    _reactant_num[i] = coupledScalar("reactants", i);
    //std::cout << (*getScalarVar("reactants", i)).name() << std::endl;
  }
}

Real
ODEReactionBase::computeProduct()
{
  _reactant_product = 1;
  for (unsigned int r = 0; r < _num_reactants; r++)
  {
    _reactant_product *= (*_reactants[r])[_i];
  }
  return _stoichiometric_coeff * _reactant_product;
}

Real
ODEReactionBase::computeProductJacobian()
{
  Real power, gas_mult;
  power = 0.0;
  // eq_u_mult = 1.0;
  gas_mult = 1.0;

  for (unsigned int r = 0; r < _num_reactants; r++)
  {
    if (_reactant_num[r] == _var.number())
    {
      power += 1;
    }
    else
    {
      gas_mult *= (*_reactants[r])[_i];
    }
  }
  return _stoichiometric_coeff * gas_mult * power * std::pow(_u[_i], power - 1);
}

Real
ODEReactionBase::computeProductOffDiagJacobian(unsigned int jvar)
{
  /*
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

  if (jvar == _w_var && !_w_eq_u)
  {
    power += 1;
    deriv_factor = _w[_i];
  }
  else
    other_factor *= _w[_i];

  if (jvar == _x_var && !_x_eq_u)
  {
    power += 1;
    deriv_factor = _x[_i];
  }
  else
    other_factor *= _x[_i];

  return -_stoichiometric_coeff * other_factor * power * std::pow(deriv_factor, power - 1);
  */
  // return 0;

  // drdv includes the derivative of the variables w.r.t. the current variable indicated by jvar
  Real power, drdv, product;
  power = 0;
  drdv = 1;
  product = 1;
  for (unsigned int r = 0; r < _num_reactants; r++)
  {
    if (_reactant_num[r] == jvar)
    {
      power += 1;
      drdv = (*_reactants[r])[_i];
    }
    else
    {
      product *= (*_reactants[r])[_i];
    }
  }
  return _stoichiometric_coeff * product * power * std::pow(drdv, power - 1);
}
