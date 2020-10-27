#pragma once

#include "ODEKernel.h"
// #include "RateCoefficientProvider.h"

class Energy2BodyScalar;

template <>
InputParameters validParams<Energy2BodyScalar>();

class Energy2BodyScalar : public ODEKernel
{
public:
  Energy2BodyScalar(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  unsigned int _v_var;
  const VariableValue & _v;
  unsigned int _w_var;
  const VariableValue & _w;
  const VariableValue & _rate_coefficient;
  const Real _threshold_energy;
};
