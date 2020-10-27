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

  /*
  unsigned int _v_var;
  const VariableValue & _v;
  const VariableValue & _rate_coefficient;
  */
};
