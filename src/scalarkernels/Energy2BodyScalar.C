#include "Energy2BodyScalar.h"

registerMooseObject("CraneApp", Energy2BodyScalar);

template <>
InputParameters
validParams<Energy2BodyScalar>()
{
  InputParameters params = validParams<ODEKernel>();
  params.addRequiredCoupledVar("v", "Coupled variable 1.");
  params.addRequiredCoupledVar("w", "Coupled variable 2.");
  params.addCoupledVar("rate_coefficient", 0, "Coupled reaction coefficient (if equation-based).");
  params.addRequiredParam<Real>("threshold_energy", "The energy change of this reaction."); 
  return params;
}

Energy2BodyScalar::Energy2BodyScalar(const InputParameters & parameters)
  : ODEKernel(parameters),
    _v_var(isCoupledScalar("v") ? coupledScalar("v") : -1),
    _v(coupledScalarValue("v")),
    _w_var(isCoupledScalar("w") ? coupledScalar("w") : -1),
    _w(coupledScalarValue("w")),
    _rate_coefficient(coupledScalarValue("rate_coefficient")),
    _threshold_energy(getParam<Real>("threshold_energy"))
{
}

Real
Energy2BodyScalar::computeQpResidual()
{
  return -_rate_coefficient[_i] * _v[_i] * _w[_i] * _threshold_energy;
}

Real
Energy2BodyScalar::computeQpJacobian()
{
  return 0;
}

Real
Energy2BodyScalar::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real power, deriv_factor, other_factor;
  power = 0;
  other_factor = 1;
  deriv_factor = 1;

  if (jvar == _v_var)
  {
    power += 1;
    deriv_factor = _v[_i];
  }
  else
    other_factor *= _v[_i];

  if (jvar == _w_var)
  {
    power += 1;
    deriv_factor = _w[_i];
  }
  else
    other_factor *= _w[_i];

  return -_rate_coefficient[_i] * other_factor * power *
         std::pow(deriv_factor, power - 1) * _threshold_energy;
}
