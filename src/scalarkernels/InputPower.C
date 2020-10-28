#include "InputPower.h"

registerMooseObject("CraneApp", InputPower);

template <>
InputParameters
validParams<InputPower>()
{
  InputParameters params = validParams<ODEKernel>();
  params.addParam<Real>("value", 1, "The input power density.");
  params.addCoupledVar(
      "power", 1, "Coupled auxiliary variable representing power (if applicable).");
  params.addClassDescription(
      "This class is intended to be used to add input power density to the electron energy "
      "equation (and optionally the gas temperature). Returns a power density of value * power, "
      "where 'value' is a constant power density and 'power' is an Auxiliary Variable. Both "
      "default to 1.");
  return params;
}

InputPower::InputPower(const InputParameters & parameters)
  : ODEKernel(parameters), _value(getParam<Real>("value")), _power(coupledScalarValue("power"))
{
}

Real
InputPower::computeQpResidual()
{
  return -_value * _power[_i] / 1.602e-19;
}

Real
InputPower::computeQpJacobian()
{
  return 0;
}

Real
InputPower::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0;
}
