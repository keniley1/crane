/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#pragma once

#include "ReactionLogBase.h"

class ConstantReactionLog : public ReactionLogBase
{
public:
  static InputParameters validParams();
  ConstantReactionLog(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  Real _r_units;

  const ADMaterialProperty<Real> & _rate_coefficient;

  Real _stoichiometric_value;
};
