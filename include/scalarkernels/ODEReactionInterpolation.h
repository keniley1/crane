#pragma once

#include "ODEReactionBase.h"
#include "LinearInterpolation.h"

//class ODEReactionInterpolation;
class ODEReactionInterpolation : public ODEReactionBase
{
public:
  static InputParameters validParams();

  ODEReactionInterpolation(const InputParameters & parameters);

  //virtual ~ODEReactionInterpolation() {}

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  void interpolateData();

  LinearInterpolation _rate_coefficient;
  const VariableValue & _sampler_var;
  Real _sampler_const;
  std::string _sampling_format;
  bool _use_time;
  bool _use_log;
  Real _scale_factor;
};
