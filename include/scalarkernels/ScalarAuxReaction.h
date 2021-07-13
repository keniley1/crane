#pragma once

#include "ScalarReactionBase.h"
// #include "RateCoefficientProvider.h"

class ScalarAuxReaction;

template <>
InputParameters validParams<ScalarAuxReaction>();

class ScalarAuxReaction : public ScalarReactionBase
{
public:
  static InputParameters validParams();
  ScalarAuxReaction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  Real _stoichiometric_coeff;
  const VariableValue & _rate_coefficient;
};
