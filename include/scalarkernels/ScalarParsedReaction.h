#pragma once

#include "ScalarReactionBase.h"
#include "FunctionParserUtils.h"
// #include "RateCoefficientProvider.h"

class ScalarParsedReaction;

template <>
InputParameters validParams<ScalarParsedReaction>();

class ScalarParsedReaction : public ScalarReactionBase, public FunctionParserUtils<false>
{
public:
  static InputParameters validParams();
  ScalarParsedReaction(const InputParameters & parameters);

protected:
  void updateParams();

  virtual Real rateCoefficient();
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// function expression
  std::string _function;

  /// coupled variables
  unsigned int _nargs;
  std::vector<const VariableValue *> _args;
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

  Real _stoichiometric_coeff;

private:
  /// Vector to look up the internal coupled variable index into _arg_*  through the libMesh variable number
  std::vector<unsigned int> _arg_index;
};
