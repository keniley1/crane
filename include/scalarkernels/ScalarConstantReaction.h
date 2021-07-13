#pragma once

#include "ScalarReactionBase.h"
// #include "RateCoefficientProvider.h"

class ScalarConstantReaction;

template <>
InputParameters validParams<ScalarConstantReaction>();

class ScalarConstantReaction : public ScalarReactionBase
{
public:
  static InputParameters validParams();
  ScalarConstantReaction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /*
  unsigned int _v_var;
  const VariableValue & _v;
  const VariableValue & _rate_coefficient;

  Real _stoichiometric_coeff;
  // Real _reaction_coeff;
  bool _v_eq_u;
  bool _rate_constant_equation;

  // const RateCoefficientProvider & _data;
  */

  Real _stoichiometric_coeff;
  Real _rate_coefficient;
};
