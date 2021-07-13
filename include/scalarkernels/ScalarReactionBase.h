#pragma once

#include "ODEKernel.h"
// #include "RateCoefficientProvider.h"

class ScalarReactionBase;

template <>
InputParameters validParams<ScalarReactionBase>();

class ScalarReactionBase : public ODEKernel
{
public:
  static InputParameters validParams();
  ScalarReactionBase(const InputParameters & parameters);

protected:
  /*
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  unsigned int _v_var;
  const VariableValue & _v;
  const VariableValue & _rate_coefficient;

  Real _stoichiometric_coeff;
  // Real _reaction_coeff;
  bool _v_eq_u;
  bool _rate_constant_equation;

  // const RateCoefficientProvider & _data;
  */
  virtual Real computeQpResidual() = 0;
  virtual Real multiplyReactants();
  virtual Real multiplyReactantsJacobian(unsigned int ivar);
  virtual Real multiplyReactantsDerivative(unsigned int ivar);

  std::vector<Real> _coefficients;
  //const bool _townsend;

  // Reactant variables
  std::vector<const VariableValue *> _reactants;
  //std::vector<const unsigned int> _reactants_var;
  std::vector<unsigned int> _reactants_var;
  unsigned int _num_reactants;

  Real _val;
};
