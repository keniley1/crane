#include "ODEReactionFunction.h"

// MOOSE includes
#include "MooseVariableScalar.h"
#include "SystemBase.h"

#include "libmesh/fparser_ad.hh"

registerMooseObject("CraneApp", ODEReactionFunction);

defineLegacyParams(ODEReactionFunction);

InputParameters
ODEReactionFunction::validParams()
{
  // InputParameters params = validParams<ODEKernel>();
  InputParameters params = ODEReactionBase::validParams();
  params += FunctionParserUtils<false>::validParams();
  params.addClassDescription("Rate equation with a rate coefficient defined by a function.");

  params.addRequiredParam<std::string>("function", "function expression");
  params.addCoupledVar("args", "additional coupled variables");
  params.addParam<std::vector<std::string>>(
      "constant_names", "Vector of constants used in the parsed function (use this for kB etc.)");
  params.addParam<std::vector<std::string>>(
      "constant_expressions",
      "Vector of values for the constants in constant_names (can be an FParser expression)");
  params.addParam<std::vector<PostprocessorName>>(
      "postprocessors", "Vector of postprocessor names used in the function expression");
  return params;
}

ODEReactionFunction::ODEReactionFunction(const InputParameters & parameters)
  : ODEReactionBase(parameters),
    FunctionParserUtils(parameters),
    _function(getParam<std::string>("function")),
    _nargs(coupledScalarComponents("args")),
    _args(_nargs),
    _arg_names(_nargs),
    _func_dFdarg(_nargs),
    _number_of_nl_variables(_sys.nVariables()),
    _arg_index(_number_of_nl_variables, -1)
{
  // build variables argument (start with variable the kernel is operating on)
  std::string variables = _var.name();

  // add additional coupled variables
  for (unsigned int i = 0; i < _nargs; ++i)
  {
    _arg_names[i] = getScalarVar("args", i)->name();
    variables += "," + _arg_names[i];
    _args[i] = &coupledScalarValue("args", i);

    // populate number -> arg index lookup table skipping aux variables
    unsigned int number = coupledScalar("args", i);
    if (number < _number_of_nl_variables)
      _arg_index[number] = i;
  }

  // add postprocessors
  auto pp_names = getParam<std::vector<PostprocessorName>>("postprocessors");
  _pp.resize(pp_names.size());
  for (unsigned int i = 0; i < pp_names.size(); ++i)
  {
    variables += "," + pp_names[i];
    _pp[i] = &getPostprocessorValueByName(pp_names[i]);
  }

  // base function object
  _func_F = std::make_shared<SymFunction>();

  // set FParser interneal feature flags
  setParserFeatureFlags(_func_F);

  // add the constant expressions
  addFParserConstants(_func_F,
                      getParam<std::vector<std::string>>("constant_names"),
                      getParam<std::vector<std::string>>("constant_expressions"));

  // parse function
  std::cout << _func_F->Parse(_function, variables) << std::endl;
  if (_func_F->Parse(_function, variables) >= 0)
    mooseError("Invalid function\n",
               _function,
               "\nin ParsedODEKernel ",
               name(),
               ".\n",
               _func_F->ErrorMsg());

  // on-diagonal derivative
  _func_dFdu = std::make_shared<SymFunction>(*_func_F);

  if (_func_dFdu->AutoDiff(_var.name()) != -1)
    mooseError("Failed to take first derivative w.r.t. ", _var.name());

  // off-diagonal derivatives
  for (unsigned int i = 0; i < _nargs; ++i)
  {
    _func_dFdarg[i] = std::make_shared<SymFunction>(*_func_F);

    if (_func_dFdarg[i]->AutoDiff(_arg_names[i]) != -1)
      mooseError("Failed to take first derivative w.r.t. ", _arg_names[i]);
  }

  // optimize
  if (!_disable_fpoptimizer)
  {
    _func_F->Optimize();
    _func_dFdu->Optimize();
    for (unsigned int i = 0; i < _nargs; ++i)
      _func_dFdarg[i]->Optimize();
  }

  // just-in-time compile
  if (_enable_jit)
  {
    _func_F->JITCompile();
    _func_dFdu->JITCompile();
    for (unsigned int i = 0; i < _nargs; ++i)
      _func_dFdarg[i]->JITCompile();
  }

  // reserve storage for parameter passing buffer
  _func_params.resize(1 + _nargs + _pp.size());
}

void
ODEReactionFunction::updateParams()
{
  _func_params[0] = _u[_i];

  for (unsigned int j = 0; j < _nargs; ++j)
    _func_params[j + 1] = (*_args[j])[_i];
  for (unsigned int j = 0; j < _pp.size(); ++j)
    _func_params[j + 1 + _nargs] = *_pp[j];
}

Real
ODEReactionFunction::computeQpResidual()
{
  // return -_rate_coefficient * computeProduct();
  updateParams();
  return -computeProduct() * evaluate(_func_F);
}

Real
ODEReactionFunction::computeQpJacobian()
{
  updateParams();
  // By the product rule, we can split the rate coefficient jacobian term and the product jacobian
  // terms up and sum them together.
  return -(computeProduct() * evaluate(_func_dFdu) + computeProductJacobian() * evaluate(_func_F));
  // return 0;
}

Real
ODEReactionFunction::computeQpOffDiagJacobian(unsigned int jvar)
{
  updateParams();
  int i = _arg_index[jvar];
  if (i < 0)
    return -(computeProductOffDiagJacobian(jvar) * evaluate(_func_F));
  else
  {
    return -(computeProduct() * evaluate(_func_dFdarg[i]) +
             computeProductOffDiagJacobian(jvar) * evaluate(_func_F));
  }
  // return 0;
}
