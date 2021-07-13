#include "ScalarAuxReaction.h"

registerMooseObject("CraneApp", ScalarAuxReaction);

template <>
InputParameters
validParams<ScalarAuxReaction>()
{
  InputParameters params = validParams<ScalarReactionBase>();
  params.addRequiredParam<Real>("nu",
                                "The number produced or destroyed in this reaction. Positive "
                                "denotes a gain, negative denotes a loss.");
  // params.addRequiredParam<Real>("rate_coefficient",
  //                              "The constant rate coefficient of this reaction.");
  params.addCoupledVar("rate_coefficient",
                       "The coupled AuxVariable representing the rate coefficient.");
  params.addClassDescription("Base class for generic reactions. Multiplies the reactants of the "
                             "specified reaction together to be supplied to the child class.");
  return params;
}

ScalarAuxReaction::ScalarAuxReaction(const InputParameters & parameters)
  : ScalarReactionBase(parameters),
    _stoichiometric_coeff(getParam<Real>("nu")),
    _rate_coefficient(coupledScalarValue("rate_coefficient"))
{
}

Real
ScalarAuxReaction::computeQpResidual()
{
  return -_stoichiometric_coeff * _rate_coefficient[_i] * ScalarReactionBase::multiplyReactants();
}

Real
ScalarAuxReaction::computeQpJacobian()
{
  return -_stoichiometric_coeff * _rate_coefficient[_i] *
         ScalarReactionBase::multiplyReactantsDerivative(_var.number());
}

Real
ScalarAuxReaction::computeQpOffDiagJacobian(unsigned int jvar)
{
  return -_stoichiometric_coeff * _rate_coefficient[_i] *
         ScalarReactionBase::multiplyReactantsDerivative(jvar);
}
