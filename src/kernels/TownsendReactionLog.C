#include "TownsendReactionLog.h"

// MOOSE includes
#include "MooseUtils.h"
#include "MooseVariable.h"

registerADMooseObject("CraneApp", TownsendReactionLog);

InputParameters
TownsendReactionLog::validParams()
{
  InputParameters params = ReactionLogBase::validParams();

  params.addParam<std::string>(
      "rate_coefficient_name",
      "The name of the material property containing this townsend coefficient. By default the "
      "Crane will add its own name, but if added manually a user-defined name may be provided.");
  params.addRequiredCoupledVar("potential", "The potential.");
  params.addRequiredCoupledVar("electrons", "The electron density.");
  params.addRequiredParam<Real>("nu",
                                "The total stoichiometric coefficient for the variable this "
                                "reaction is acting on. This is positive or negative, depending on "
                                "whether the variable is produced or consumed by this reaction.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredParam<std::string>("reaction", "Stores the full reaction equation.");
  params.addParam<std::string>(
      "number",
      "",
      "The reaction number. Optional, just for material property naming purposes. If a single "
      "reaction has multiple different rate coefficients (frequently the case when multiple "
      "species are lumped together to simplify a reaction network), this will prevent the same "
      "material property from being declared multiple times.");
  return params;
}

TownsendReactionLog::TownsendReactionLog(const InputParameters & parameters)
  : ReactionLogBase(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _diffem(getADMaterialProperty<Real>("diffem")),
    _muem(getADMaterialProperty<Real>("muem")),
    _townsend_coefficient(getADMaterialProperty<Real>(
        isParamValid("rate_coefficient_name")
            ? getParam<std::string>("rate_coefficient_name")
            : "alpha" + getParam<std::string>("number") + "_" + getParam<std::string>("reaction"))),
    _grad_potential(adCoupledGradient("potential")),
    _em(adCoupledValue("electrons")),
    _stoichiometric_value(getParam<Real>("nu")),
    _grad_em(adCoupledGradient("electrons"))
{
}

ADReal
TownsendReactionLog::computeQpResidual()
{
  return -_test[_i][_qp] * ReactionLogBase::multiplyReactants() *
         (std::exp(_em[_qp]) * (-_muem[_qp] * -_grad_potential[_qp] * _r_units -
                                _diffem[_qp] * _grad_em[_qp] * _r_units))
             .norm() *
         _townsend_coefficient[_qp] * _stoichiometric_value;

  /*
  return -_test[_i][_qp] * _alpha[_qp] * std::exp(_target[_qp]) *
         (std::exp(_em[_qp]) * (-_muem[_qp] * -_grad_potential[_qp] * _r_units -
                                _diffem[_qp] * _grad_em[_qp] * _r_units))
             .norm() *
         _coefficient;
         */
}
