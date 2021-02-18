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

// Parrot includes
#include "FlowAndTransport.h"
#include "ParrotSolver.h"

// #include "MeshModifier.h"
// #include "FractureUserObject.h"



// #include "libmesh/quadrature.h"

// Forward declarations
class SolveDiffusion;

template <>
InputParameters validParams<SolveDiffusion>();

class SolveDiffusion : public GeneralUserObject
{
public:
  static InputParameters validParams();

  SolveDiffusion(const InputParameters & params);

  virtual void initialize ()                   override;
  virtual void execute    ()                   override {};
  virtual void finalize   ()                   override {};
  virtual void threadJoin (const UserObject &) override {};

  protected:

  	std::string  _var_name;
    std::string  _sys_name;
    std::string  _material_name;

    unsigned int _p_var;

	void init();
	void assemble();
	void setDirichlet();
	void solve();
	void write();
	void del();

  QBase const * const & _qrule;

  //Parallel::Communicator * _pp_comm;
  EquationSystems        * _equationSystemsM;
  MeshBase               * _meshBase;
	EquationSystems        * _equationSystemsP;
	LinearImplicitSystem   * _linearImplicitSystemP;
	DistributedMesh        * _mesh;

	ParrotSolver * _parrotSolver;

  std::vector< Parallel::Communicator const * > _pp_comm;

	unsigned int _dim;

	DofMap * _dof_map;

	MaterialBase     * _materialBase;
	FlowAndTransport * _flowAndTransport;

	std::vector<BoundaryName>     _dirichletNodesetNames;
  std::vector<BoundaryID>       _dirichletSidesetIds;
	std::vector<Function const *> _dirichletFunctions;
	std::vector<FunctionName>     _dirichletFunctionNames;

  std::vector<BoundaryName>     _neumannSidesetNames;
  std::vector<BoundaryID>       _neumannSidesetIds;
  std::vector<Function const *> _neumannFunctions;
  std::vector<FunctionName>     _neumannFunctionNames;

	int _solverType;

	bool const _has_exodus_file;
	bool const _has_nemesis_file;
  bool const _has_variable_name;

	std::string _exodus_filename;
	std::string _nemesis_filename;
  std::string _variable_name;

  std::vector<dof_id_type>  _dirIds;
    // std::vector<AuxVariableName> _aux_var_names;

    // AuxVariableName _aux_var_name;

    // std::vector<int> _vector_p;
    // std::vector<Real> _vector_value;
    // std::vector<boundary_id_type> _boundary_D_ids;
    // std::vector<boundary_id_type> _boundary_N_ids;
    // std::vector<Real> _value_N_bc;
    // std::vector<Real> _value_D_bc;

    // void AssembleDiffusionOP(EquationSystems & _es, const std::string & system_name);

    // Real ComputeMaterialProprties(const Elem *elem);

    // int solve(EquationSystems & _es);

    // void set_solution(EquationSystems & _es);

};

