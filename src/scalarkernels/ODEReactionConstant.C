#include "ODEReactionConstant.h"

registerMooseObject("CraneApp", ODEReactionConstant);

defineLegacyParams(ODEReactionConstant);

InputParameters
ODEReactionConstant::validParams()
{
  // InputParameters params = validParams<ODEKernel>();
  InputParameters params = ODEReactionBase::validParams();
  params.addRequiredParam<Real>("rate_coefficient", "Rate coefficient for the specified reaction.");
  return params;
}

ODEReactionConstant::ODEReactionConstant(const InputParameters & parameters)
  : ODEReactionBase(parameters), _rate_coefficient(getParam<Real>("rate_coefficient"))
{
}

Real
ODEReactionConstant::computeQpResidual()
{
  return -_rate_coefficient * computeProduct();
}

Real
ODEReactionConstant::computeQpJacobian()
{
  return -_rate_coefficient * computeProductJacobian();
}

Real
ODEReactionConstant::computeQpOffDiagJacobian(unsigned int jvar)
{
  return -_rate_coefficient * computeProductOffDiagJacobian(jvar);
}
