#pragma once

#include "ODEKernel.h"
// #include "RateCoefficientProvider.h"

class InputPower;

template <>
InputParameters validParams<InputPower>();

class InputPower : public ODEKernel
{
public:
  InputPower(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const Real _value;
  const VariableValue & _power;
};
