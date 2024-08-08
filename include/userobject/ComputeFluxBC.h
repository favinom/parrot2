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

class ComputeFluxBC;

template <>
InputParameters validParams<ComputeFluxBC>();

class ComputeFluxBC : public GeneralUserObject
{
public:
  static InputParameters validParams();

  ComputeFluxBC(const InputParameters & params);

  virtual void initialize ()                   override;
  virtual void execute    ()                   override {};
  virtual void finalize   ()                   override {};
  virtual void threadJoin (const UserObject &) override {};

  protected:

	void init();
	void assemble();
	void setDirichlet();
	void solve();
	void write();
	void del();

  std::string  _material_name;
  std::string  _var_name;
  std::string  _sys_name;
  std::string  _mat_domain_name_0;
  // std::string  _mat_binterp_name_0;
  // std::string  _mat_binterp_name_1;
  std::string  _vec_flux_name_0;
  // std::string  _vec_flux_name_1;
  std::string  _vec_rhs_name_0;
  // std::string  _vec_rhs_name_1;
  // std::string  _vec_boundary_name_0;
  // std::string  _vec_boundary_name_1;
  std::string  _vec_dir;

  QBase const * const & _qrule;

  // int _block_0_id;
  // int _block_1_id;

  EquationSystems       & _equationSystemsM;
  MeshBase        const * _meshBase;
	LinearImplicitSystem  * _linearImplicitSystemM;


  unsigned int _v_var;

  MaterialBase     * _materialBase;
  FlowAndTransport * _flowAndTransport;

	unsigned int _dim;

	DofMap * _dof_map;

  std::vector<BoundaryName>     _dirichletNodesetNames;
  std::vector<BoundaryID>       _dirichletSidesetIds;
	std::vector<Function const *> _dirichletFunctions;
	std::vector<FunctionName>     _dirichletFunctionNames;

  std::vector<BoundaryName>     _neumannSidesetNames;
  std::vector<BoundaryID>       _neumannSidesetIds;
  std::vector<Function const *> _neumannFunctions;
  std::vector<FunctionName>     _neumannFunctionNames;

	bool const _has_exodus_file;
	// bool const _has_nemesis_file;
  // bool const _has_variable_name_1;

  // bool const _has_variable_name_2;

	std::string _exodus_filename;
	// std::string _nemesis_filename;
  // std::string _variable_name_1;

  // std::string _variable_name_2;

  std::vector<dof_id_type>  _dirIds;

  

};

