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
  //std::cout << -_rate_coefficient * computeProduct() << std::endl;
  return -_rate_coefficient * computeProduct();
  // std::cout << _rate_coefficient << ", " << (*_reactants[0])[_i] << ", " << (*_reactants[1])[_i]
  // << ", " << _stoichiometric_coeff << std::endl; return -_rate_coefficient * (*_reactants[0])[_i]
  // * (*_reactants[1])[_i] * _stoichiometric_coeff;
  //std::cout << -_rate_coefficient * _stoichiometric_coeff * (*_reactants[0])[_i] * (*_reactants[0])[_i] << std::endl;
  //return -_rate_coefficient * _stoichiometric_coeff * (*_reactants[0])[_i] * (*_reactants[0])[_i];
  //std::cout << -_rate_coefficient * _stoichiometric_coeff * _u[_i] * _u[_i] << std::endl;
  //return -_rate_coefficient  * _stoichiometric_coeff * _u[_i] * _u[_i];
}

Real
ODEReactionConstant::computeQpJacobian()
{
  /*
  Real power, eq_u_mult, gas_mult;
  power = 0.0;
  eq_u_mult = 1.0;
  gas_mult = 1.0;

  if (isCoupledScalar("v") && _v_eq_u)
  {
    power += 1.0;
    eq_u_mult = _v[_i];
  }
  else
    gas_mult *= _v[_i];

  if (isCoupledScalar("w") && _w_eq_u)
  {
    power += 1.0;
    eq_u_mult = _w[_i];
  }
  else
    gas_mult *= _w[_i];

  if (isCoupledScalar("x") && _x_eq_u)
  {
    power += 1.0;
    eq_u_mult = _x[_i];
  }
  else
    gas_mult *= _x[_i];

  return -_stoichiometric_coeff * _rate_coefficient[_i] * gas_mult * power *
         std::pow(eq_u_mult, power - 1);
         */
  // return 0;

  return -_rate_coefficient * computeProductJacobian();
  //return 0;
}

Real
ODEReactionConstant::computeQpOffDiagJacobian(unsigned int jvar)
{
  /*
  Real power, deriv_factor, other_factor;
  power = 0;
  other_factor = 1;
  deriv_factor = 1;

  if (jvar == _v_var && !_v_eq_u)
  {
    power += 1;
    deriv_factor = _v[_i];
  }
  else
    other_factor *= _v[_i];

  if (jvar == _w_var && !_w_eq_u)
  {
    power += 1;
    deriv_factor = _w[_i];
  }
  else
    other_factor *= _w[_i];

  if (jvar == _x_var && !_x_eq_u)
  {
    power += 1;
    deriv_factor = _x[_i];
  }
  else
    other_factor *= _x[_i];

  return -_stoichiometric_coeff * _rate_coefficient[_i] * other_factor * power *
         std::pow(deriv_factor, power - 1);
         */
  return -_rate_coefficient * computeProductOffDiagJacobian(jvar);
  //return 0;
}
