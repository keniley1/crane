#include "ReactionLogBase.h"

// MOOSE includes
#include "MooseUtils.h"
#include "MooseVariable.h"

//registerADMooseObject("CraneApp", ReactionLogBase);

InputParameters
ReactionLogBase::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addParam<std::string>(
      "electron_name", "em", "The name of the electrons. Zapdos automatically defaults to `em`.");
  params.addParam<bool>("townsend",
                        false,
                        "Whether or not this is reaction based on a Townsend coefficient. If true, "
                        "electrons should not be included in the reactant list.");
  params.addRequiredCoupledVar("reactants",
                               "All of the reactants involved in this reaction. If this is a "
                               "townsend reaction, electrons are excluded.");
  params.addRequiredParam<std::vector<Real>>(
      "coefficients",
      "The number of each reactant that participates in this reaction. For example, in the "
      "ionization reaction, e + Ar -> 2e + Ar+, the coefficient of both e and Ar is 1 because only "
      "one of each appears on the reactant side of the equation.");
  params.addRequiredParam<std::string>("reaction", "Stores the full reaction equation.");
  params.addParam<std::string>(
      "number",
      "",
      "The reaction number. Optional, just for material property naming purposes. If a single "
      "reaction has multiple different rate coefficients (frequently the case when multiple "
      "species are lumped together to simplify a reaction network), this will prevent the same "
      "material property from being declared multiple times.");
  params.addClassDescription("Base class for generic reactions. Multiplies the reactants of the "
                             "specified reaction together to be supplied to the child class.");
  return params;
}

ReactionLogBase::ReactionLogBase(const InputParameters & parameters)
  : ADKernel(parameters),
    _coefficients(getParam<std::vector<Real>>("coefficients")),
    _townsend(getParam<bool>("townsend"))
{
  int dval = -1;
  _num_reactants = coupledComponents("reactants");

  if (_townsend && _num_reactants > 2)
    mooseError("Only two body Townsend reactions are supported. Please make sure that the reaction "
               "is written correctly, and if so, set townsend = false for this reaction.");

  _reactants.resize(_num_reactants);

  for (unsigned int i = 0; i < _num_reactants; ++i)
  {
    // If this is a townsend reaction, check to make sure electrons are
    // not included in the reactant list. Those are multiplied in later.
    // If included, it is ignored and that element is erased from
    // the reactant vector.
    // Foolproof! (Hopefully...)
    if (_townsend && getVar("reactants", i)->name() == getParam<std::string>("electron_name"))
    {
      dval = i;
      continue;
    }

    _reactants[i] = &adCoupledValue("reactants", i);
  }

  // Delete the element corresponding to electrons from reactant list if this
  // is labeled as a Townsend reaction.
  if (dval >= 0)
  {
    if (_coefficients.size() == _reactants.size())
      _coefficients.erase(_coefficients.begin() + dval);

    _reactants.erase(_reactants.begin() + dval);
    _num_reactants += -1;
  }

  if (_num_reactants != _coefficients.size())
    mooseError("The number of reactants and coefficients is not equal!");
}

ADReal
ReactionLogBase::multiplyReactants()
{
  _val = 0;
  for (unsigned int i = 0; i < _num_reactants; ++i)
  {
    _val += (*_reactants[i])[_qp] * _coefficients[i];
  }

  return std::exp(_val);
}
// ReactionLogBase::computeQpResidual()
