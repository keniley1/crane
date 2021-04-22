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

#include "ADKernel.h"

class ReactionLogBase : public ADKernel
{
public:
  static InputParameters validParams();
  ReactionLogBase(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  std::vector<Real> _coefficients;
  const bool _townsend;

  // Reactant variables
  std::vector<const ADVariableValue *> _reactants;
  unsigned int _num_reactants;

  ADReal _val;
};
