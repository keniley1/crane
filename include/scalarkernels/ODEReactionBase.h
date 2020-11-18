#pragma once

#include "ODEKernel.h"
// #include "RateCoefficientProvider.h"

class ODEReactionBase : public ODEKernel
{
public:
  static InputParameters validParams();

  ODEReactionBase(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() = 0;
  virtual Real computeQpJacobian() = 0;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) = 0;

  virtual Real computeProduct();
  virtual Real computeProductJacobian();
  virtual Real computeProductOffDiagJacobian(unsigned int jvar);

  unsigned int _num_reactants;
  unsigned int _reactant_var;
  std::vector<VariableValue *> _reactants;
  std::vector<int> _reactant_num;
  // const VariableValue & _rate_coefficient;

  Real _stoichiometric_coeff;

  Real _reactant_product;
  // bool _rate_constant_equation;

  // const RateCoefficientProvider & _data;
};
