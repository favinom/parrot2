//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PermeabilityDiffusion.h"

registerMooseObject("parrot2App", PermeabilityDiffusion);

defineLegacyParams(PermeabilityDiffusion);

InputParameters
PermeabilityDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  return params;
}

PermeabilityDiffusion::PermeabilityDiffusion(const InputParameters & parameters) :
Kernel(parameters),
_Kscalar(getMaterialProperty<Real>("conductivityProperty"))
{}

Real
PermeabilityDiffusion::computeQpResidual()
{

  return _Kscalar[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
PermeabilityDiffusion::computeQpJacobian()
{
  return _Kscalar[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
