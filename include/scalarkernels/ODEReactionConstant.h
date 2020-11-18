#pragma once

#include "ODEReactionBase.h"

//class ODEReactionConstant;
class ODEReactionConstant : public ODEReactionBase
{
public:
  static InputParameters validParams();

  ODEReactionConstant(const InputParameters & parameters);

  //virtual ~ODEReactionConstant() {}

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const Real _rate_coefficient;
};
