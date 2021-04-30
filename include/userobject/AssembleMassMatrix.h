//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "StoreOperators.h"
// Parrot includes
#include "FlowAndTransport.h"
#include "libmesh/linear_implicit_system.h"

// Forward declarations
class AssembleMassMatrix;

template <>
InputParameters validParams<AssembleMassMatrix>();

class AssembleMassMatrix : public GeneralUserObject
{
public:
 AssembleMassMatrix(const InputParameters & params);

  virtual void initialize() override;
  virtual void execute() override ;
  virtual void finalize() override {};
  virtual void threadJoin(const UserObject & ) override {};
  static void getRow(PetscMatrix<Number> & matrix, int const & row, std::vector<Real> & values, std::vector<int> & columns)
{
    Mat const & mat=matrix.mat();
    PetscInt ncol;
    PetscInt const *col;
    PetscScalar const *val;
    MatGetRow(mat,row,&ncol,&col,&val);
    values.resize(ncol);
    columns.resize(ncol);
    for (int i=0; i<ncol; ++i)
    {
        values[i] =val[i];
        columns[i]=col[i];
    }
    MatRestoreRow(mat,row,&ncol,&col,&val);
}
  protected:

    // std::vector<int> _vector_p;
    // std::vector<Real> _vector_value;
    std::string  _material_name;

    std::vector< Parallel::Communicator const * > _pp_comm;


    unsigned int _test_var;

    UserObjectName const _userObjectName;
    StoreOperators * _storeOperatorsUO;
    EquationSystems        * _equationSystemsM;
    MeshBase               * _meshBase;
    EquationSystems        * _equationSystemsT;
    LinearImplicitSystem   * _linearImplicitSystemT;
    DistributedMesh        * _mesh;

    // std::shared_ptr<SparseMatrix<Number>> _interpolator;
    // std::shared_ptr<PetscMatrix<Number>> _mass_matrix;
    // std::shared_ptr<PetscMatrix<Number>> _lump_mass_matrix;
    // std::shared_ptr<PetscMatrix<Number>> _poro_mass_matrix;
    // std::shared_ptr<PetscMatrix<Number>> _poro_lump_mass_matrix;
    // std::shared_ptr<PetscMatrix<Number>> _hanging_interpolator;
    std::shared_ptr<PetscVector<Number>> _hanging_vec;

    void assemble_mass_matrix();
    // Real ComputeMaterialProprties(const Elem *elem);


    MaterialBase     * _materialBase;
    FlowAndTransport * _flowAndTransport;

    bool const _constrainMatrices;
    bool const _code_dof_map;

    bool _hasMeshModifier;
    std::string _meshModifierName;

    QBase const * const & _qrule;


    std::shared_ptr<PetscVector<Number>> _bc_vec;

    std::shared_ptr<PetscVector<Number>> _value_bc_vec;
 
    std::string _dc_var;

    std::vector<int> _dc_boundary_id;

    std::vector<std::vector<int> > _dc_variables_id;

    std::vector<Real> _value_D_bc;

    void find_boundary(std::vector<int> &_dc_boundary_id);

    void determine_dc_bnd_var_id(const std::vector<std::string> & BC_var);

    std::vector<std::string> split_string(const std::string & s, char delim);
};

