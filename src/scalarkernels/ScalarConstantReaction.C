#include "ScalarConstantReaction.h"

registerMooseObject("CraneApp", ScalarConstantReaction);

template <>
InputParameters
validParams<ScalarConstantReaction>()
{
  InputParameters params = validParams<ScalarReactionBase>();
  params.addRequiredParam<Real>(
      "nu",
      "The number produced or destroyed in this reaction. Positive denotes a gain, negative denotes a loss.");
  params.addRequiredParam<Real>(
      "rate_coefficient",
      "The constant rate coefficient of this reaction.");
  params.addClassDescription("Base class for generic reactions. Multiplies the reactants of the "
                             "specified reaction together to be supplied to the child class.");
  return params;
}

ScalarConstantReaction::ScalarConstantReaction(const InputParameters & parameters)
  : ScalarReactionBase(parameters),
  /*
    _v_var(isCoupledScalar("v") ? coupledScalar("v") : -1),
    _v(coupledScalarValue("v")),
    _rate_coefficient(coupledScalarValue("rate_coefficient")),
    _stoichiometric_coeff(getParam<Real>("coefficient")),
    _v_eq_u(getParam<bool>("v_eq_u")),
    _rate_constant_equation(getParam<bool>("rate_constant_equation"))
    */
    _stoichiometric_coeff(getParam<Real>("nu")),
    _rate_coefficient(getParam<Real>("rate_coefficient"))
{
}


Real
ScalarConstantReaction::computeQpResidual()
{
  return -_stoichiometric_coeff * _rate_coefficient * ScalarReactionBase::multiplyReactants();
}

Real
ScalarConstantReaction::computeQpJacobian()
{
  return -_stoichiometric_coeff * _rate_coefficient * ScalarReactionBase::multiplyReactantsDerivative(_var.number());
}

Real
ScalarConstantReaction::computeQpOffDiagJacobian(unsigned int jvar)
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

  return -_stoichiometric_coeff * _rate_coefficient[_i] * other_factor * power *
         std::pow(deriv_factor, power - 1);
         */
  return -_stoichiometric_coeff * _rate_coefficient * ScalarReactionBase::multiplyReactantsDerivative(jvar);
}
