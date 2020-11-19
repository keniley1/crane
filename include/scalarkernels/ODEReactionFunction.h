#pragma once

#include "ODEReactionBase.h"
#include "FunctionParserUtils.h"

// class ODEReactionFunction;
class ODEReactionFunction : public ODEReactionBase, public FunctionParserUtils<false>
{
public:
  static InputParameters validParams();

  ODEReactionFunction(const InputParameters & parameters);

  // virtual ~ODEReactionFunction() {}

protected:
  void updateParams();
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// function expression
  std::string _function;

  /// coupled variables
  unsigned int _nargs;
  std::vector<VariableValue *> _args;
  std::vector<std::string> _arg_names;

  /// function parser object for the residual and on-diagonal Jacobian
  SymFunctionPtr _func_F;
  SymFunctionPtr _func_dFdu;

  /// function parser objects for the Jacobian
  std::vector<SymFunctionPtr> _func_dFdarg;

  /// number of non-linear variables in the problem
  const unsigned int _number_of_nl_variables;

  /// coupled postprocessors
  std::vector<const PostprocessorValue *> _pp;

  usingFunctionParserUtilsMembers(false);


private:
  /// Vector to look up the internal coupled variable index into _arg_*  through the libMesh variable number
  std::vector<unsigned int> _arg_index;
};