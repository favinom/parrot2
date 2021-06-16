//* This file is part of the MOOSE framework

//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

// Forward Declarations
class ElementIntegralConcentrationPostprocessor;

template <>
InputParameters validParams<ElementIntegralConcentrationPostprocessor>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class ElementIntegralConcentrationPostprocessor : public ElementIntegralPostprocessor,
                                             public MooseVariableInterface<Real>
{
public:
  ElementIntegralConcentrationPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Holds the solution at current quadrature points
 
   MeshGeneratorName       _meshGeneratorName;
   bool              const _hasMeshGenerator;
   unsigned int  _regionId;
   const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
//  const VariableGradient & _grad_u;
};


