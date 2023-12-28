//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssignMonomialValue.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseMeshUtils.h"

// LibMesh includes

#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"


#include "helpers.h"



registerMooseObject("parrot2App", AssignMonomialValue);

defineLegacyParams(AssignMonomialValue);

InputParameters AssignMonomialValue::validParams()
{
  InputParameters params = GeneralUserObject::validParams();\

  params.addParam<std::string>("variable_name", "the variable name of the output");


  return params;
}

AssignMonomialValue::AssignMonomialValue(const InputParameters & parameters) :
GeneralUserObject(parameters),
_has_variable_name( isParamValid("variable_name") )


{


  if (_has_variable_name)
    _variable_name=getParam<std::string>("variable_name");
}




void AssignMonomialValue::initialize()
{
   

    MeshBase & mesh = _fe_problem.mesh();
    
    MooseVariable & _var = _fe_problem.getStandardVariable(0,  _variable_name);
    
    System & _sys = _var.sys().system();

    NumericVector<Number> * _solution = _sys.solution.get();
    
    // PetscInt       rstart,rend;
    // PetscScalar    tmp_sol;
    // VecGetOwnershipRange(utopia::raw_type(sol),&rstart,&rend);
    
    for (const auto & elem : mesh.active_local_element_ptr_range())
    for (unsigned int comp = 0; comp < elem->n_comp(_sys.number(), _var.number()); comp++)
    {
      const dof_id_type _index = elem->dof_number(_sys.number(), _var.number(), comp);
      _solution->set(_index, elem->id());
    }

  _solution->close();
  _sys.update();
    
}


