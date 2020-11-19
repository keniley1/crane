#pragma once

#include "AddVariableAction.h"
#include "Action.h"
#include "ChemicalReactionsBase.h"

class AddScalarReactions;

template <>
InputParameters validParams<AddScalarReactions>();

// class ChemicalReactions : public AddVariableAction
class AddScalarReactions : public ChemicalReactionsBase
{
public:
  AddScalarReactions(InputParameters params);
  virtual void addSuperelasticRateCoefficient(const unsigned & reaction_num);
  virtual void addEEDFEnergy(const unsigned & reaction_num, const std::string & kernel_name);
  virtual void addFunctionKernel(const unsigned & reaction_num,
                                 const unsigned & species_num,
                                 const std::string & kernel_name,
                                 const bool & energy_kernel);
  virtual void addConstantKernel(const unsigned & reaction_num,
                                 const unsigned & species_num,
                                 const std::string & kernel_name,
                                 const bool & energy_kernel);
  virtual std::string
  getKernelName(const unsigned & num_reactants, const bool & energy_kernel, const bool & is_aux);

  const std::string _interpolation_type;
  // AddScalarReactions(InputParameters params) : ChemicalReactionsBase(params) {};

  virtual void act();

  // protected:
  std::vector<std::string> _aux_scalar_var_name;

  std::string _log_append;
};
