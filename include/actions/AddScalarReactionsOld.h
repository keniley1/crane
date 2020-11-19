#pragma once

#include "AddVariableAction.h"
#include "Action.h"
#include "ChemicalReactionsBase.h"

class AddScalarReactionsOld;

template <>
InputParameters validParams<AddScalarReactionsOld>();

// class ChemicalReactions : public AddVariableAction
class AddScalarReactionsOld : public ChemicalReactionsBase
{
public:
  AddScalarReactionsOld(InputParameters params);
  const std::string _interpolation_type;
  // AddScalarReactionsOld(InputParameters params) : ChemicalReactionsBase(params) {};

  virtual void act();

//protected:
  std::vector<std::string> _aux_scalar_var_name;
};
