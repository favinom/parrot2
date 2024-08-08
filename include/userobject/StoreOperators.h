///* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
//#include "UserObject.h"

#include "GeneralUserObject.h"
#include "libmesh/petsc_matrix.h"

#include "libmesh/petsc_vector.h"

// Forward declarations
class StoreOperators;

template <>
InputParameters validParams<StoreOperators>();


class StoreOperators : public GeneralUserObject
{
public:
  StoreOperators(const InputParameters & params);
  virtual void initialize() override {};
  virtual void execute() override
  {
//    _mass_matrix->print_matlab("massC.txt");
//    _interpolator->print_matlab("int.txt");
//    exit(1);
  };
  virtual void finalize() override   {};
  virtual void threadJoin(const UserObject & ) override {};
 
  std::shared_ptr<PetscMatrix<Number>> &
  Interpolator()
  {
    return _interpolator;
  };

  std::shared_ptr<PetscMatrix<Number>> &
  H_Interpolator()
  {
    return _hanging_interpolator;
  };

  std::shared_ptr<PetscMatrix<Number>> &
  MassMatrix()
  {
    return _mass_matrix;
  };


  std::shared_ptr<PetscMatrix<Number>> &
  PoroMassMatrix()
  {
    return _poro_mass_matrix;
  };

  std::shared_ptr<PetscMatrix<Number>> &
  LumpMassMatrix()
  {
    return _lump_mass_matrix;
  };

  std::shared_ptr<PetscMatrix<Number>> &
  PoroLumpMassMatrix()
  {
    return _poro_lump_mass_matrix;
  };

 std::shared_ptr<PetscMatrix<Number>> &
  JacMatrix()
  {
    return _jac_matrix;
  };

  PetscMatrix<Number> * &
  StabMatrix()
  {
    return _stab_matrix;
  };

  std::shared_ptr<PetscVector<Number>> &
  BcVec()
  {
    return _bc_vec;
  };

  std::shared_ptr<PetscVector<Number>> &
  ValueBcVec()
  {
    return _value_bc_vec;
  };

  std::shared_ptr<PetscVector<Number>> &
  SolVec()
  {
    return _solution_vec;
  };

  std::shared_ptr<PetscVector<Number>> &
  HangVec()
  {
    return _hanging_vec;
  };




protected:
  std::shared_ptr<PetscMatrix<Number>> _interpolator;
  std::shared_ptr<PetscMatrix<Number>> _hanging_interpolator;
  std::shared_ptr<PetscVector<Number>> _hanging_vec;
  std::shared_ptr<PetscMatrix<Number>> _mass_matrix;
  std::shared_ptr<PetscMatrix<Number>> _poro_mass_matrix;
  std::shared_ptr<PetscMatrix<Number>> _lump_mass_matrix;
  std::shared_ptr<PetscMatrix<Number>>  _poro_lump_mass_matrix;
  std::shared_ptr<PetscMatrix<Number>>  _jac_matrix;
  std::shared_ptr<PetscVector<Number>>  _bc_vec;
  std::shared_ptr<PetscVector<Number>>  _value_bc_vec;
  std::shared_ptr<PetscVector<Number>>  _solution_vec;
  PetscMatrix<Number> * _stab_matrix;

};




