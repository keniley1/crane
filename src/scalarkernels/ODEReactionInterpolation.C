#include "ODEReactionInterpolation.h"

registerMooseObject("CraneApp", ODEReactionInterpolation);

defineLegacyParams(ODEReactionInterpolation);

InputParameters
ODEReactionInterpolation::validParams()
{
  InputParameters params = ODEReactionBase::validParams();
  params.addCoupledVar("sampler", 0, "The variable with which the data will be sampled.");
  params.addParam<bool>("use_time", false, "Whether or not to sample with time.");
  params.addParam<bool>(
      "use_log", false, "Whether or not to return the natural logarithm of the sampled data.");
  params.addParam<Real>("scale_factor",
                        1.0,
                        "Multiplies the sampled output by a given factor. Convert from m^3 to "
                        "cm^3, for example.(Optional)");
  params.addParam<Real>("const_sampler", 0, "The value with which the data will be sampled.");
  params.addParam<std::string>(
      "property_file", "", "The file containing interpolation tables for material properties.");
  params.addParam<std::string>(
      "file_location", "", "The name of the file that stores the reaction rate tables.");
  params.addParam<std::string>("sampling_format",
                               "reduced_field",
                               "The format that the rate constant files are in. Options: "
                               "reduced_field and electron_energy.");
  return params;
}

ODEReactionInterpolation::ODEReactionInterpolation(const InputParameters & parameters)
  : ODEReactionBase(parameters),
    _sampler_var(coupledScalarValue("sampler")),
    _sampler_const(getParam<Real>("const_sampler")),
    _sampling_format(getParam<std::string>("sampling_format")),
    _use_time(getParam<bool>("use_time")),
    _use_log(getParam<bool>("use_log")),
    _scale_factor(getParam<Real>("scale_factor"))
{
  // interpolateData is wrapped into a function because in the future, this table might be
  // updated dynamically.
  // This will probably be done through a UserObject that this class will optionally accept and
  // read data from.
  interpolateData();
}

void
ODEReactionInterpolation::interpolateData()
{
  std::vector<Real> x_val;
  std::vector<Real> y_val;
  std::string file_name =
      getParam<std::string>("file_location") + "/" + getParam<std::string>("property_file");
  MooseUtils::checkFileReadable(file_name);
  const char * charPath = file_name.c_str();
  std::ifstream myfile(charPath);
  Real value;

  if (myfile.is_open())
  {
    while (myfile >> value)
    {
      x_val.push_back(value);
      myfile >> value;
      y_val.push_back(value);
    }
    myfile.close();
  }
  else
    mooseError("Unable to open file");

  _rate_coefficient.setData(x_val, y_val);
}

Real
ODEReactionInterpolation::computeQpResidual()
{
  return -_rate_coefficient.sample(_sampler_var[_i]) * computeProduct();
}

Real
ODEReactionInterpolation::computeQpJacobian()
{
  // return -_rate_coefficient * computeProductJacobian();
  return -_rate_coefficient.sample(_sampler_var[_i]) * computeProductJacobian();
  // return 0;
}

Real
ODEReactionInterpolation::computeQpOffDiagJacobian(unsigned int jvar)
{
  // return -_rate_coefficient * computeProductOffDiagJacobian(jvar);
  return -_rate_coefficient.sample(_sampler_var[_i]) * computeProductOffDiagJacobian(jvar);
  // return 0;
}
