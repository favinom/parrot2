//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "GeneralUserObject.h"
#include "Function.h"

// LibMesh includes
#include "libmesh/linear_implicit_system.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/mesh.h"
#include "libmesh/distributed_mesh.h"


class AssignMonomialValue;

template <>
InputParameters validParams<AssignMonomialValue>();

class AssignMonomialValue: public GeneralUserObject
{
public:
  static InputParameters validParams();

  AssignMonomialValue(const InputParameters & params);

  virtual void initialize ()                   override;
  virtual void execute    ()                   override {};
  virtual void finalize   ()                   override {};
  virtual void threadJoin (const UserObject &) override {};

  protected:

  	std::string  _var_name;
    std::string  _sys_name;
    std::string  _material_name;

    
  bool const _has_variable_name;
    std::string _variable_name;

};

