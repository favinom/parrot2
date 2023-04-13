//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeFluxUtopia.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseMeshUtils.h"
#include "MultiApp.h"

// LibMesh includes
#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "helpers.h"

#include "utopia_MeshTransferOperator.hpp"

inline void convert_vec(libMesh::NumericVector<libMesh::Number> &lm_vec,  utopia::UVector &utopia_vec) {
    using namespace libMesh;
    Vec p_vec = cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&lm_vec)->vec();
    utopia_vec.wrap(p_vec);
}

inline void convert_mat(libMesh::SparseMatrix<libMesh::Number> &lm_mat, utopia::USparseMatrix &utopia_mat) {
    using namespace libMesh;

    Mat p_mat = cast_ptr<libMesh::PetscMatrix<libMesh::Number> *>(&lm_mat)->mat();
    utopia::convert(p_mat, utopia_mat);
    }

typedef utopia::USparseMatrix SparseMatT;


registerMooseObject("parrot2App", ComputeFluxUtopia);

defineLegacyParams(ComputeFluxUtopia);

InputParameters ComputeFluxUtopia::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<int>("block_0_id", "block_0_id");

  params.addRequiredParam<int>("block_1_id", "block_1_id");

  params.addParam<int>("block_2_id", "block_2_id");
      
  params.addRequiredParam< std::vector<BoundaryName> >("dirichlet_nodeset_names_M","The name of the nodeset to create");

  params.addRequiredParam< std::vector<FunctionName> >("dirichlet_function_names_M", "The function.");

    params.addRequiredParam< std::vector<BoundaryName> >("dirichlet_nodeset_names_F","The name of the nodeset to create");

  params.addRequiredParam< std::vector<FunctionName> >("dirichlet_function_names_F", "The function.");

  params.addParam< std::vector<BoundaryName> >("neumann_sideset_names",std::vector<BoundaryName>(),"The ids of the nodeset to create");

  params.addParam< std::vector<FunctionName> >("neumann_function_names",std::vector<FunctionName>(), "The function.");

  params.addParam<std::string>("exodus_file", "the file name of the output");

  params.addParam<std::string>("nemesis_file", "the file name of the output");

  params.addParam<std::string>("variable_name_1", "the variable name of the output");

  params.addParam<std::string>("variable_name_2", "the variable name of the output");

  params.addRequiredParam<std::vector<int>>("material_block_id","material_block_id");

  params.addRequiredParam<std::vector<Real>>("material_value_p","material_value_p");

  params.addRequiredParam<MultiAppName>("multi_app", "The MultiApp's name in your input file!");
  return params;

}

ComputeFluxUtopia::ComputeFluxUtopia(const InputParameters & parameters) :
GeneralUserObject(parameters),
_qrule(_assembly.qRule()),
_block_0_id( getParam< int >( "block_0_id" ) ),
_block_1_id( getParam< int >( "block_1_id" ) ),
_equationSystemsM( _fe_problem.es() ),
_dirichletNodesetNames_M ( getParam<std::vector<BoundaryName> >("dirichlet_nodeset_names_M" ) ),
_dirichletFunctionNames_M( getParam<std::vector<FunctionName> >("dirichlet_function_names_M") ),
_dirichletNodesetNames_F ( getParam<std::vector<BoundaryName> >("dirichlet_nodeset_names_F" ) ),
_dirichletFunctionNames_F( getParam<std::vector<FunctionName> >("dirichlet_function_names_F") ),
_neumannSidesetNames   ( getParam<std::vector<BoundaryName> >("neumann_sideset_names"   ) ),
_neumannFunctionNames  ( getParam<std::vector<FunctionName> >("neumann_function_names"  ) ),
_has_exodus_file(  isParamValid("exodus_file")  ),
_has_variable_name_1( isParamValid("variable_name_1") ),
_has_variable_name_2( isParamValid("variable_name_2") ),
_has_block_2_id(isParamValid( "block_2_id" )),
_vector_p(getParam<std::vector<int>>("material_block_id")),
_vector_value(getParam<std::vector<Real>>("material_value_p")),
_multiapp_name(getParam<MultiAppName>("multi_app"))

{
  _sys_name_M           ="Flux_system_M";
  _var_name_M           ="Dummy_variable_M";
  _mat_domain_name_M0   ="matrix_domain_M0";
  _mat_domain_name_M1   ="matrix_domain_M1";
  _mat_domain_name_M2   ="matrix_domain_M2";
  _vec_flux_name_M0     ="vector_flux_M0";
  _vec_flux_name_M1     ="vector_flux_M1";
  _vec_flux_name_M2     ="vector_flux_M2";
  _vec_rhs_name_M0      ="vector_rhs_M0";
  _vec_rhs_name_M1      ="vector_rhs_M1";
  _vec_rhs_name_M2      ="vector_rhs_M2";
  _vec_boundary_name_M0 ="vec_boundary_M0";
  _vec_boundary_name_M1 ="vec_boundary_M1";
  _vec_boundary_name_M2 ="vec_boundary_M2";
  _vec_dir_M            ="vector_dirichlet_M";



  _sys_name_F           ="Flux_system_F";
  _var_name_F           ="Dummy_variable_F";
  _mat_domain_name_F0   ="matrix_domain_F0";
  _mat_domain_name_F1   ="matrix_domain_F1";
  _mat_domain_name_F2   ="matrix_domain_F2";
  _vec_flux_name_F0     ="vector_flux_F0";
  _vec_flux_name_F1     ="vector_flux_F1";
  _vec_flux_name_F2     ="vector_flux_F2";
  _vec_rhs_name_F0      ="vector_rhs_F0";
  _vec_rhs_name_F1      ="vector_rhs_F1";
  _vec_rhs_name_F2      ="vector_rhs_F2";
  _vec_boundary_name_F0 ="vec_boundary_F0";
  _vec_boundary_name_F1 ="vec_boundary_F1";
  _vec_boundary_name_F2 ="vec_boundary_F2";
  _vec_dir_F            ="vector_dirichlet_F";

  if (_has_exodus_file)
    _exodus_filename=getParam<std::string>("exodus_file");
  
  if (_has_variable_name_1){
    _variable_name_1=getParam<std::string>("variable_name_1");
  }

  if (_has_variable_name_2){
    _variable_name_2=getParam<std::string>("variable_name_2");
  }

  if ( _dirichletNodesetNames_M.size() != _dirichletFunctionNames_M.size() )
  {
    mooseError("_dirichletNodesetNames_M.size() != _dirichletFunctionNames_M.size()");
  }

  if ( _dirichletNodesetNames_F.size() != _dirichletFunctionNames_F.size() )
  {
    mooseError("_dirichletNodesetNames_F.size() != _dirichletFunctionNames_F.size()");
  }


  if ( _neumannSidesetNames.size() != _neumannFunctionNames.size() )
  {
    mooseError("_neumannSidesetNames.size() != _neumannFunctionNames.size()");
  }

  if (_has_block_2_id){
    _block_2_id =  getParam< int >( "block_2_id" );
  }
}

void ComputeFluxUtopia::init()
{
  std::cout<<"ComputeFluxUtopia::init() start\n";

  _meshBase = &_equationSystemsM.get_mesh();

  MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
    
  FEProblemBase & F_problem = _multi_app.appProblemBase(0);

  _meshBaseF  = &F_problem.mesh().getMesh();

  _linearImplicitSystemM = &_equationSystemsM.add_system<LinearImplicitSystem> (_sys_name_M.c_str());
  
  _v_varM = _linearImplicitSystemM[0].add_variable (_var_name_M.c_str(), FIRST);
      
  _linearImplicitSystemM[0].add_matrix(_mat_domain_name_M0.c_str());
  _linearImplicitSystemM[0].add_matrix(_mat_domain_name_M1.c_str());
  _linearImplicitSystemM[0].add_matrix(_mat_domain_name_M2.c_str());

  _linearImplicitSystemM[0].add_vector(_vec_flux_name_M0.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_flux_name_M1.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_flux_name_M2.c_str());

  _linearImplicitSystemM[0].add_vector(_vec_rhs_name_M0.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_rhs_name_M1.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_rhs_name_M2.c_str());

  _linearImplicitSystemM[0].add_vector(_vec_boundary_name_M0.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_boundary_name_M1.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_boundary_name_M2.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_dir_M.c_str());


  MooseVariable & _f_var_p = F_problem.getStandardVariable(0,"pressure_F");

  EquationSystems &_equationSystemsF = _f_var_p.sys().system().get_equation_systems();


  _linearImplicitSystemF =&_equationSystemsF.add_system<LinearImplicitSystem> (_sys_name_F.c_str());

  _v_varF = _linearImplicitSystemF[0].add_variable (_var_name_F.c_str(), FIRST);

  _linearImplicitSystemF[0].add_matrix(_mat_domain_name_F0.c_str());
  _linearImplicitSystemF[0].add_matrix(_mat_domain_name_F1.c_str());
  _linearImplicitSystemF[0].add_matrix(_mat_domain_name_F2.c_str());

  _linearImplicitSystemF[0].add_vector(_vec_flux_name_F0.c_str());
  _linearImplicitSystemF[0].add_vector(_vec_flux_name_F1.c_str());
  _linearImplicitSystemF[0].add_vector(_vec_flux_name_F2.c_str());

  _linearImplicitSystemF[0].add_vector(_vec_rhs_name_F0.c_str());
  _linearImplicitSystemF[0].add_vector(_vec_rhs_name_F1.c_str());
  _linearImplicitSystemF[0].add_vector(_vec_rhs_name_F2.c_str());

  _linearImplicitSystemF[0].add_vector(_vec_boundary_name_F0.c_str());
  _linearImplicitSystemF[0].add_vector(_vec_boundary_name_F1.c_str());
  _linearImplicitSystemF[0].add_vector(_vec_boundary_name_F2.c_str());
  _linearImplicitSystemF[0].add_vector(_vec_dir_F.c_str());

  std::cout<<"Before reinit, it may take a while...\n";
  _equationSystemsM.reinit();
  _equationSystemsF.reinit();

  std::cout<<"Done\n";
  _equationSystemsM.print_info();
  _equationSystemsF.print_info();

  _dim_M = _meshBase[0].mesh_dimension();
  _dim_F = _meshBaseF[0].mesh_dimension();

  _dof_map = &_linearImplicitSystemM[0].get_dof_map();

  _dof_mapF = &_linearImplicitSystemF[0].get_dof_map();

  for (int i=0; i<_dirichletFunctionNames_M.size(); ++i)
  {
    Function const & func=getFunctionByName( _dirichletFunctionNames_M.at(i) );
    _dirichletFunctions_M.push_back(&func);
  }


  for (int i=0; i<_dirichletFunctionNames_F.size(); ++i)
  {
    Function const & func_F=getFunctionByName( _dirichletFunctionNames_F.at(i) );
    _dirichletFunctions_F.push_back(&func_F);
  }

  for (int i=0; i<_neumannFunctionNames.size(); ++i)
  {
    Function const & func=getFunctionByName( _neumannFunctionNames.at(i) );
    _neumannFunctions.push_back(&func);
  }

  _dirichletSidesetIds_M = MooseMeshUtils::getBoundaryIDs(_meshBase[0], _dirichletNodesetNames_M, true);
  _dirichletSidesetIds_F = MooseMeshUtils::getBoundaryIDs(_meshBaseF[0], _dirichletNodesetNames_F, true);
  _neumannSidesetIds   = MooseMeshUtils::getBoundaryIDs(_meshBase[0], _neumannSidesetNames,   true);

  std::cout<<"ComputeFluxUtopia::init() stop\n";
}

void ComputeFluxUtopia::setDirichletMatrix()
{
  std::cout<<"ComputeFluxUtopia::setDirichletMatrix() start\n";

  NumericVector<Number> & dM = _linearImplicitSystemM[0].get_vector(_vec_dir_M.c_str());
  dM=1.0;


  BoundaryInfo const & boundary_info = _meshBase[0].get_boundary_info();
  for ( auto & node : as_range(_meshBase->active_nodes_begin(), _meshBase->active_nodes_end()) )
  {
    for (int i=0; i<_dirichletSidesetIds_M.size(); ++i)
    {
      if ( boundary_info.has_boundary_id(node,_dirichletSidesetIds_M.at(i)) )
      {
        std::vector<dof_id_type> di;
        _dof_map[0].dof_indices (node,di, 0);
        if (di.size()!=1)
        {
          mooseError("di.size()!=1");
        }
        dof_id_type index=di.at(0);
        _dirIds.push_back(index);

        dM.set(index,0.0);
        break;
      }
    }
  }

  dM.close();


  std::cout<<"ComputeFluxUtopia::setDirichletMatrix() stop\n";
}


void ComputeFluxUtopia::setDirichletFracture()
{
  std::cout<<"ComputeFluxUtopia::setDirichletFracture() start\n";

  NumericVector<Number> & dF = _linearImplicitSystemF[0].get_vector(_vec_dir_F.c_str());
 
  dF=1.0;

  BoundaryInfo const & boundary_info = _meshBaseF[0].get_boundary_info();

  for ( auto & nodeF : as_range(_meshBaseF->active_nodes_begin(), _meshBaseF->active_nodes_end()) )
  {
    for (int i=0; i<_dirichletSidesetIds_F.size(); ++i)
    {
     // std::cout<<"dirichletBoundar is "<<_dirichletSidesetIds_F.at(i)<<"and size is "<< _dirichletSidesetIds_F.size()<<std::endl;

      if ( boundary_info.has_boundary_id(nodeF,_dirichletSidesetIds_F.at(i)) )
      {
        std::vector<dof_id_type> di;
        _dof_mapF[0].dof_indices (nodeF,di, 0);
        if (di.size()!=1)
        {
          mooseError("di.size()!=1");
        }
        dof_id_type index=di.at(0);
        _dirIds.push_back(index);
        //std::cout<<"index"<<index<<" and "<<_dirichletSidesetIds_F.at(i)<<std::endl;
        dF.set(index,0.0);
        break;
      }
    }
  }

  dF.close();


  std::cout<<"ComputeFluxUtopia::setDirichletFracture() stop\n";
}

void ComputeFluxUtopia::solve()
{
  std::cout<<"ComputeFluxUtopia::solveMatrix() start\n";

  SparseMatrix<Number> & AM0 = _linearImplicitSystemM[0].get_matrix(_mat_domain_name_M0.c_str());
  SparseMatrix<Number> & AM1 = _linearImplicitSystemM[0].get_matrix(_mat_domain_name_M1.c_str());
  SparseMatrix<Number> & AM2 = _linearImplicitSystemM[0].get_matrix(_mat_domain_name_M2.c_str());

  NumericVector<Number> & fM0 = _linearImplicitSystemM[0].get_vector(_vec_flux_name_M0.c_str());
  NumericVector<Number> & fM1 = _linearImplicitSystemM[0].get_vector(_vec_flux_name_M1.c_str());
  NumericVector<Number> & fM2 = _linearImplicitSystemM[0].get_vector(_vec_flux_name_M2.c_str());
  
  NumericVector<Number> & bM0 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_M0.c_str());
  NumericVector<Number> & bM1 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_M1.c_str());
  NumericVector<Number> & bM2 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_M2.c_str());
  
  NumericVector<Number> & dM  = _linearImplicitSystemM[0].get_vector(_vec_dir_M.c_str());

  NumericVector<Number> & iM0 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_M0.c_str());
  NumericVector<Number> & iM1 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_M1.c_str());
  NumericVector<Number> & iM2 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_M2.c_str());


 {


    MooseVariableFEBase  & _sol_aux = _fe_problem.getVariable(0, "aux_sol_M", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
    MooseVariableFEBase  & _sol_main = _fe_problem.getVariable(0, "pressure_M", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

    
    System & main_sys = _sol_main.sys().system();
    System & aux_sys = _sol_aux.sys().system();
    NumericVector<Number> * aux_solution = aux_sys.solution.get();
    NumericVector<Number> * main_solution = main_sys.solution.get();

    {

      for (const auto & node : _fe_problem.es().get_mesh().local_node_ptr_range())
      {

        for (unsigned int comp = 0; comp < node->n_comp(main_sys.number(), _sol_main.number()); comp++)
          
          {
            
            const dof_id_type proj_index = node->dof_number(main_sys.number(), _sol_main.number(), comp);

            const dof_id_type to_index_1 = node->dof_number(aux_sys.number(), _sol_aux.number(), comp);

            main_solution->set(proj_index, ((*aux_solution)(to_index_1)));
            }
          }
        }

        main_solution->close();
        main_sys.update();
      }


    
      LinearImplicitSystem & linearImplicitSystem = _equationSystemsM.get_system<LinearImplicitSystem> ("nl0");
      std::unique_ptr<NumericVector<Number>> & solutionPointer(linearImplicitSystem.solution);
      NumericVector<Number> * solM = solutionPointer.get();

      AM0.vector_mult(fM0,solM[0]);
      fM0*=-1.0;
      fM0+=bM0;
      fM0*=dM;
      fM0.close();


      AM1.vector_mult(fM1,solM[0]);
      fM1*=-1.0;
      fM1+=bM1;
      fM1*=dM;
      fM1.close();


      // utopia::UVector rhs_m0;
      // utopia::USparseMatrix A_m0;

      // utopia::convert(vec_m0, rhs_m0);
      // utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*mat_m0).mat(), A_m0);
      // libMesh::PetscMatrix<libMesh::Number> *mat_m0 = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(&A0);
      // Vec vec_m0 = static_cast<libMesh::PetscVector<libMesh::Number> *>(&f0)->vec();

      if (_has_block_2_id){
        AM2.vector_mult(fM2,solM[0]);
        fM2*=-1.0;
        fM2+=bM2;
        fM2*=dM;
        fM2.close();
      }

      std::cout<<fM0.sum()<<std::endl;
      std::cout<<fM1.sum()<<std::endl;

      if(_has_block_2_id) std::cout<<fM2.sum()<<std::endl;

      std::cout<<iM0.sum()<<std::endl;
      std::cout<<iM1.sum()<<std::endl;

      if(_has_block_2_id) std::cout<<iM2.sum()<<std::endl;

      dM=fM0;
      dM*=iM0;

      bM0.close();

      dM=fM1;
      dM*=iM1;
      bM1.close();

      if(_has_block_2_id){
      dM=fM2;
      dM*=iM2;
      bM2.close();
    }

  std::cout<<"ComputeFluxUtopia::solveMatrix() stop\n";


  std::cout<<"ComputeFluxUtopia::solveFracture() start\n";

  SparseMatrix<Number> & AF0 = _linearImplicitSystemF[0].get_matrix(_mat_domain_name_F0.c_str());
  SparseMatrix<Number> & AF1 = _linearImplicitSystemF[0].get_matrix(_mat_domain_name_F1.c_str());
  SparseMatrix<Number> & AF2 = _linearImplicitSystemF[0].get_matrix(_mat_domain_name_F2.c_str());

  NumericVector<Number> & fF0 = _linearImplicitSystemF[0].get_vector(_vec_flux_name_F0.c_str());
  NumericVector<Number> & fF1 = _linearImplicitSystemF[0].get_vector(_vec_flux_name_F1.c_str());
  NumericVector<Number> & fF2 = _linearImplicitSystemF[0].get_vector(_vec_flux_name_F2.c_str());
  
  NumericVector<Number> & bF0 = _linearImplicitSystemF[0].get_vector(_vec_rhs_name_F0.c_str());
  NumericVector<Number> & bF1 = _linearImplicitSystemF[0].get_vector(_vec_rhs_name_F1.c_str());
  NumericVector<Number> & bF2 = _linearImplicitSystemF[0].get_vector(_vec_rhs_name_F2.c_str());
  
  NumericVector<Number> & dF  = _linearImplicitSystemF[0].get_vector(_vec_dir_F.c_str());

  NumericVector<Number> & iF0 = _linearImplicitSystemF[0].get_vector(_vec_boundary_name_F0.c_str());
  NumericVector<Number> & iF1 = _linearImplicitSystemF[0].get_vector(_vec_boundary_name_F1.c_str());
  NumericVector<Number> & iF2 = _linearImplicitSystemF[0].get_vector(_vec_boundary_name_F2.c_str());

  MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
  FEProblemBase & F_problem = _multi_app.appProblemBase(0);


 {
    MooseVariableFEBase  & _sol_auxF  = F_problem.getVariable(0, "aux_sol_F", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
    MooseVariableFEBase  & _sol_mainF = F_problem.getVariable(0, "pressure_F", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

    
    System & main_sysF = _sol_mainF.sys().system();
    System & aux_sysF = _sol_auxF.sys().system();
    NumericVector<Number> * aux_solutionF = aux_sysF.solution.get();
    NumericVector<Number> * main_solutionF = main_sysF.solution.get();

  {

      for (const auto & node : F_problem.es().get_mesh().local_node_ptr_range())
      {

        for (unsigned int comp = 0; comp < node->n_comp(main_sysF.number(), _sol_mainF.number()); comp++)
          
        {
            
            const dof_id_type proj_indexF= node->dof_number(main_sysF.number(), _sol_mainF.number(), comp);

            const dof_id_type to_index_1F = node->dof_number(aux_sysF.number(), _sol_auxF.number(), comp);

            main_solutionF->set(proj_indexF, ((*aux_solutionF)(to_index_1F)));
        }
      }
    }

      main_solutionF->close();
      main_sysF.update();
  }

  MooseVariable & _f_var = F_problem.getStandardVariable(0,"pressure_F");
  EquationSystems &_equationSystemsF = _f_var.sys().system().get_equation_systems();
  LinearImplicitSystem & linearImplicitSystemF = _equationSystemsF.get_system<LinearImplicitSystem> ("nl0");
  std::unique_ptr<NumericVector<Number>> & solutionPointerF(linearImplicitSystemF.solution);
  NumericVector<Number> * solF = solutionPointerF.get();

  AF0.vector_mult(fF0,solF[0]);
  fF0*=-1.0;
  fF0+=bF0;
  fF0*=dF;
  fF0.close();


  AF1.vector_mult(fF1,solF[0]);
  fF1*=-1.0;
  fF1+=bF1;
  fF1*=dF;
  fF1.close();

  if (_has_block_2_id){
    AF2.vector_mult(fF2,solF[0]);
    fF2*=-1.0;
    fF2+=bF2;
    fF2*=dF;
    fF2.close();
  }

  std::cout<<fF0.sum()<<std::endl;
  std::cout<<fF1.sum()<<std::endl;

  if(_has_block_2_id) std::cout<<fM2.sum()<<std::endl;

  std::cout<<iF0.sum()<<std::endl;
  std::cout<<iF1.sum()<<std::endl;

  if(_has_block_2_id) std::cout<<iM2.sum()<<std::endl;


  std::cout<<"ComputeFluxUtopia::solveFracture() stop\n";

  std::cout<<"Total Flux side 0 is "<<fF0.sum() + fM0.sum()<<std::endl;
  std::cout<<"Total Flux side 1 is "<<fF1.sum() + fM1.sum()<<std::endl;

  utopia::TransferOptions opts;
       
  opts.from_var_num     =  _v_varM; //master_system.variable_number("var_1_s"); //_master_aux.number();

  opts.to_var_num       =  _v_varF; //slave_system.variable_number("var_1_f"); //_slave_aux.number();

  //opts.from_var_num     = 0; // _var_1_from; //_master_aux.number();

  //opts.to_var_num       = 0; //_slave_aux.number();

  opts.n_var            = 1;

  opts.tags             = {};

  auto cs = std::make_shared<utopia::MeshTransferOperator>(
     utopia::make_ref(_meshBase[0]),
     utopia::make_ref(*_dof_map),
     utopia::make_ref(_meshBaseF[0]),
     utopia::make_ref(*_dof_mapF),
     opts);

  cs->initialize("PSEUDO_L2_PROJECTION");

  // assert(!cs->has_operator());

  // auto &matrices = cs->matrices();
  // assert(matrices.size() == 2);

  // auto _T_ptr = matrices[0];

  // utopia::USparseMatrix A_M0, A_F0;

  // utopia::UVector rhs_M0, rhs_F0, sol_M;

  // //NonlinearSystemBase & _nl_from = _fe_problem.getNonlinearSystemBase();
    
  // //utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_from.RHS()), rhs_m);

  // convert_vec(fM0, rhs_M0);
  // convert_mat(AM0, A_M0);

  // convert_vec(fF0, rhs_F0);
  // convert_mat(AF0, A_F0);

  // convert_vec(*solM, sol_M);

  // auto _T_t   = utopia::transpose(*_T_ptr);
  // auto flux_0 = (A_M0 + _T_t * A_F0 * (*_T_ptr)) * sol_M - rhs_M0 - _T_t * rhs_F0;

}

void ComputeFluxUtopia::write()
{
      std::cout<<"ComputeFluxUtopia::write() start\n";

      if(_has_exodus_file)
      {
        ExodusII_IO (_equationSystemsM.get_mesh()).write_equation_systems(_exodus_filename.c_str(), _equationSystemsM); //output_filename.c_str()
      }

      std::cout<<"ComputeFluxUtopia::write() stop\n";
}

void ComputeFluxUtopia::del()
{
      std::cout<<"ComputeFluxUtopia::del() start\n";

      // //_equationSystemsP[0].delete_system(_sys_name.c_str());
      // _equationSystemsP[0].clear();

      std::cout<<"ComputeFluxUtopia::del() stop\n";  
}

void ComputeFluxUtopia::initialize()
{

      init();
      assembleMatrix();
      assembleFracture();
      setDirichletMatrix();
      setDirichletFracture();
      solve();
      write();
      
      //_equationSystemsM.reinit();
      
      if(_has_variable_name_1 && _has_variable_name_2)

      {
        //copyVariableFl(_equationSystemsM, _var_name, _equationSystemsM, _variable_name);

        NumericVector<Number> & b0 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_M0.c_str());
        NumericVector<Number> & b1 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_M1.c_str());

        MooseVariableFEBase  & _flux_var_1 = _fe_problem.getVariable(0, _variable_name_1, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
        MooseVariableFEBase  & _flux_var_2 = _fe_problem.getVariable(0, _variable_name_2, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
        MooseVariableFEBase  & sol_var = _fe_problem.getVariable(0, "pressure_M", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

        // solution of the original system
        
        System & main_sys = sol_var.sys().system();
        System & aux_sys = _flux_var_1.sys().system();
        NumericVector<Number> * aux_solution = aux_sys.solution.get();

        {

          for (const auto & node : _fe_problem.es().get_mesh().local_node_ptr_range())
          {

            for (unsigned int comp = 0; comp < node->n_comp(main_sys.number(), sol_var.number()); comp++)
              
              {
                
                const dof_id_type proj_index = node->dof_number(main_sys.number(), sol_var.number(), comp);

                //const dof_id_type to_index = node->dof_number(main_sys.number(), main_var.number(), comp);

                const dof_id_type to_index_1 = node->dof_number(aux_sys.number(), _flux_var_1.number(), comp);

                const dof_id_type to_index_2 = node->dof_number(aux_sys.number(), _flux_var_2.number(), comp);

                //main_solution->set(to_index, sol(proj_index));

                aux_solution->set(to_index_1, b0(proj_index));
                aux_solution->set(to_index_2, b1(proj_index));
                }
              }
            }
    
            //main_solution->close();
            aux_solution->close();
            //main_sys.update();
            aux_sys.update();

            // _console << "END CopyFlux"  << std::endl;
          }
      //del();

      std::cout<<"END initialize\n";
}

void ComputeFluxUtopia::assembleMatrix()
{
  std::cout << "ComputeFluxUtopia::assemble Matrix start\n";

  SparseMatrix <Number> & AM0=_linearImplicitSystemM[0].get_matrix( _mat_domain_name_M0.c_str() );
  SparseMatrix <Number> & AM1=_linearImplicitSystemM[0].get_matrix( _mat_domain_name_M1.c_str() );
  SparseMatrix <Number> & AM2=_linearImplicitSystemM[0].get_matrix( _mat_domain_name_M2.c_str() );

  NumericVector<Number> & bM0=_linearImplicitSystemM[0].get_vector( _vec_rhs_name_M0.c_str() );
  NumericVector<Number> & bM1=_linearImplicitSystemM[0].get_vector( _vec_rhs_name_M1.c_str() );
  NumericVector<Number> & bM2=_linearImplicitSystemM[0].get_vector( _vec_rhs_name_M2.c_str() );

  NumericVector<Number> & iM0 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_M0.c_str());
  NumericVector<Number> & iM1 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_M1.c_str());
  NumericVector<Number> & iM2 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_M2.c_str());

  AM0.zero();
  AM1.zero();
  AM2.zero();

  bM0.zero();
  bM1.zero();
  bM2.zero();

  iM0.zero();
  iM1.zero();
  iM2.zero();

  FEType fe_type = _dof_map[0].variable_type(0);

  std::unique_ptr<FEBase> fe (FEBase::build(_dim_M, fe_type));
  std::unique_ptr<QBase> qrule( QBase::build (_qrule->type(),_dim_M,_qrule->get_order()));
  fe->attach_quadrature_rule (qrule.get());

  std::unique_ptr<FEBase> fe_face (FEBase::build(_dim_M, fe_type));
  std::unique_ptr<QBase> qface(fe_type.default_quadrature_rule(_dim_M-1));
  fe_face->attach_quadrature_rule (qface.get());

  std::vector<Real>                       const & JxW = fe->get_JxW();
  std::vector<std::vector<Real> >         const & phi = fe->get_phi();
  std::vector<std::vector<RealGradient> > const & dphi = fe->get_dphi();
  std::vector<Point>                      const & q_points = fe->get_xyz();

  std::vector<Real>                       const & JxW_face      = fe_face->get_JxW();
  std::vector<std::vector<Real> >         const & phi_face      = fe_face->get_phi();
  //std::vector<std::vector<RealGradient> > const & dphi_face     = fe_face->get_dphi();
  std::vector<Point>                      const & q_points_face = fe_face->get_xyz();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_Z;
  std::vector<dof_id_type> dof_indices_B;
  // std::vector<dof_id_type> dof_indices_I;

  DenseMatrix<Number> ke;
  DenseMatrix<Number> ke_Z;
  // DenseMatrix<Number> ke_I;
  DenseVector<Number> re;

  MeshBase::const_element_iterator           el = _meshBase[0].active_local_elements_begin();
  MeshBase::const_element_iterator const end_el = _meshBase[0].active_local_elements_end();

  std::vector<Number> permeability;

  int count_0=0;

  int count_1=0;

  int count_2=0;

  for ( ; el != end_el; ++el)
  {
    Elem const * elem = *el;
    fe->reinit (elem);

    int elemId=elem->subdomain_id();

    _dof_map[0].dof_indices(elem, dof_indices);
    _dof_map[0].dof_indices(elem, dof_indices_Z);
    unsigned int const n_dofs = cast_int<unsigned int>(dof_indices.size());

    ke.resize (n_dofs , n_dofs);
    ke.zero();

    ke_Z.resize (n_dofs , n_dofs);
    ke_Z.zero();

    re.resize (n_dofs );
    re.zero();

  Real localPermeability=ComputeMaterialProprties(elem);
    
  permeability.assign( qrule->n_points() ,localPermeability );
    
  for (unsigned int i=0; i<phi.size(); i++)
    {
      for (unsigned int j=0; j<phi.size(); j++)
      {
        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        {
          ke(i,j) +=  JxW[qp] * ( dphi[j][qp] * ( permeability.at(qp) *  dphi[i][qp] ) );
        }
      }
    }

  for (auto side : elem->side_index_range())
  {
      if (elem->neighbor_ptr(side) == nullptr)
      {
        fe_face->reinit(elem, side);
        for(int k =0; k < _neumannSidesetIds.size(); k++)
        {
          if ( _meshBase[0].get_boundary_info().has_boundary_id (elem, side, _neumannSidesetIds.at(k)) )
          {
            for (unsigned int qp=0; qp<qface->n_points(); qp++)
              for (unsigned int i=0; i != n_dofs; i++){

               // std::cout<<" _neumannSidesetIds.at("<<k<<")"<<"is "<< _neumannSidesetIds.at(k)<<std::endl;
                
                re(i) += JxW_face[qp] * _neumannFunctions.at(k)[0].value(0.0,q_points_face[qp]) * phi_face[i][qp];
              }
          }
        }
      }
  }

  _dof_map[0].constrain_element_matrix_and_vector (ke, re, dof_indices,false);
  _dof_map[0].constrain_element_matrix (ke_Z, dof_indices_Z,false);

  if (ke_Z.m()>n_dofs)
  {
    for (int ii=0; ii<ke_Z.m(); ++ii)
      if (ke_Z(ii,ii)>0.5)
      {
        ke(ii,ii)=0.0;
      }
  }

  for (auto side : elem->side_index_range())
  {
    Elem const * otherElem=elem->neighbor_ptr(side);
    if (otherElem != nullptr)
    {
      int otherID=otherElem->subdomain_id();
      if (elemId!=otherID)
      {
        //std::cout<<"close to the boundary\n";
        //fe_face->reinit(elem, side);

        std::unique_ptr<const Elem> test=elem->build_side_ptr (side);
        Elem const * ciao=test.get();

        _dof_map[0].dof_indices(ciao, dof_indices_B);

        std::unique_ptr<FEBase> fel  (FEBase::build(_dim_M-1, fe_type));
        std::unique_ptr<QBase> qrulel( QBase::build (_qrule->type(),_dim_M-1,_qrule->get_order()));
        fel->attach_quadrature_rule (qrulel.get());

        std::vector<Real>                       const & JxWl = fel->get_JxW();
        std::vector<std::vector<Real> >         const & phil = fel->get_phi();
        std::vector<Point>                      const & q_pointsl = fel->get_xyz();

        fel->reinit(ciao);
        
        DenseMatrix<Number> kel;
        DenseVector<Number> rel;
        kel.resize(phil.size(),phil.size());
        rel.resize(phil.size());
        kel.zero();
        rel.zero();

        for (int i=0; i<phil.size(); ++i)
          for (int j=0; j<phil.size(); ++j)
            for (int qp=0; qp<phil.at(0).size(); ++qp)
              rel(i)+=JxWl.at(qp)*phil.at(i).at(qp)*phil.at(j).at(qp);

        _dof_map[0].constrain_element_matrix_and_vector (kel,rel, dof_indices_B,true);

        for (int i=0; i<kel.m();++i)
        {
          Real & aa=kel(i,i);
          if (aa<0.5)
            aa=1.0;
          else
          {
            aa=0.0;
            for (int j=0; j<kel.n(); ++j)
            {
              Real & bb=kel(i,j);
              if (bb<-1e-10)
              {
                bb=-1.0*bb;
              }
            }
          }
        }

        if (elemId==_block_0_id)
        {
          iM0.add_vector (rel  , dof_indices_B  );
        }
        else if(elemId==_block_1_id)
        {
          iM1.add_vector (rel  , dof_indices_B  );
        }
        else if (elemId==_block_2_id)
        {
          iM2.add_vector (rel  , dof_indices_B  );
        }
        else
        {
          mooseError("elemId does not belong to any block");
        }



      }
        
      }
  }






  if (elemId ==_block_0_id)
  {
    count_0=count_0+1;
    AM0.add_matrix (ke  , dof_indices  );
    bM0.add_vector (re  , dof_indices  );
  }
  else if (elemId ==_block_1_id)
  {
    count_1 = count_1+1;
    AM1.add_matrix (ke  , dof_indices  );
    bM1.add_vector (re  , dof_indices  );
  }
  else if (_has_block_2_id == true && elemId ==_block_2_id)
  {
    count_2 = count_2+1;
    AM2.add_matrix (ke  , dof_indices  );
    bM2.add_vector (re  , dof_indices  );
  }
  else
  {
    mooseError("elemId does not belong to any block");
  }

  }

  AM0.close();
  AM1.close();
  AM2.close();
  bM0.close();
  bM1.close();
  bM2.close();
  iM0.close();
  iM1.close();
  iM2.close();

  for (int i=iM0.first_local_index(); i<iM0.last_local_index(); ++i)
  {
    if (iM0(i)>1e-10)
    {
      iM0.set( i,1./iM0(i) );
    }
    if (iM1(i)>1e-10)
    {
      iM1.set(i,1./iM1(i));
    }
    if (iM2(i)>1e-10)
    {
      iM2.set(i,1./iM2(i));
    }

  }

  iM0.close();
  iM1.close();
  iM2.close();

  std::cout << "ComputeFluxUtopia::assemble Matrix stop\n";
}



void ComputeFluxUtopia::assembleFracture()
{
  std::cout << "ComputeFluxUtopia::assemble Fracture start\n";

  SparseMatrix <Number> & AF0=_linearImplicitSystemF[0].get_matrix( _mat_domain_name_F0.c_str() );
  SparseMatrix <Number> & AF1=_linearImplicitSystemF[0].get_matrix( _mat_domain_name_F1.c_str() );
  SparseMatrix <Number> & AF2=_linearImplicitSystemF[0].get_matrix( _mat_domain_name_F2.c_str() );

  NumericVector<Number> & bF0=_linearImplicitSystemF[0].get_vector( _vec_rhs_name_F0.c_str() );
  NumericVector<Number> & bF1=_linearImplicitSystemF[0].get_vector( _vec_rhs_name_F1.c_str() );
  NumericVector<Number> & bF2=_linearImplicitSystemF[0].get_vector( _vec_rhs_name_F2.c_str() );

  NumericVector<Number> & iF0 = _linearImplicitSystemF[0].get_vector(_vec_boundary_name_F0.c_str());
  NumericVector<Number> & iF1 = _linearImplicitSystemF[0].get_vector(_vec_boundary_name_F1.c_str());
  NumericVector<Number> & iF2 = _linearImplicitSystemF[0].get_vector(_vec_boundary_name_F2.c_str());

  AF0.zero();
  AF1.zero();
 
  AF2.zero();

  bF0.zero();
  bF1.zero();
  
  bF2.zero();

  iF0.zero();
  iF1.zero();
  
  iF2.zero();

  //FEType fe_type = _dof_mapF[0].variable_type(0);



  libMesh::FEType fe_type = _dof_mapF[0].variable_type(_v_varF);

  //fe_type.order = static_cast<libMesh::Order>(fe_type.order);

  std::unique_ptr<FEBase> fe_F (FEBase::build(_dim_F, fe_type));
  
  std::unique_ptr<QBase> qrule( QBase::build (QGAUSS,_dim_F,THIRD));
  
  fe_F->attach_quadrature_rule (qrule.get());



  std::unique_ptr<FEBase> fe_face_F (FEBase::build(_dim_F, fe_type));
  std::unique_ptr<QBase> qface(fe_type.default_quadrature_rule(_dim_F-1));
  fe_face_F->attach_quadrature_rule (qface.get());

  std::vector<Real>                       const & JxW = fe_F->get_JxW();
  std::vector<std::vector<Real> >         const & phi = fe_F->get_phi();
  std::vector<std::vector<RealGradient> > const & dphi = fe_F->get_dphi();
  std::vector<Point>                      const & q_points = fe_F->get_xyz();

  std::vector<Real>                       const & JxW_face      = fe_face_F->get_JxW();
  std::vector<std::vector<Real> >         const & phi_face      = fe_face_F->get_phi();
  std::vector<Point>                      const & q_points_face = fe_face_F->get_xyz();



  std::cout<<"_dim_F==>"<<_dim_F<<std::endl;

  std::cout<<THIRD<<std::endl;

  std::cout<<QGAUSS<<std::endl;

  std::vector<dof_id_type> dof_indices_F;
  std::vector<dof_id_type> dof_indices_FZ;
  std::vector<dof_id_type> dof_indices_FB;


  DenseMatrix<Number> ke;
  DenseMatrix<Number> ke_Z;
  DenseVector<Number> re;

  MeshBase::const_element_iterator           el = _meshBaseF[0].active_local_elements_begin();
  MeshBase::const_element_iterator const end_el = _meshBaseF[0].active_local_elements_end();

  std::vector<Number> permeability;

  int count_0=0;

  int count_1=0;

  int count_2=0;

  for ( ; el != end_el; ++el)
  {
    Elem const * elem = *el;
    fe_F->reinit (elem);

    int elemId=elem->subdomain_id();

    _dof_mapF[0].dof_indices(elem, dof_indices_F);
    _dof_mapF[0].dof_indices(elem, dof_indices_FZ);

    unsigned int const n_dofs = cast_int<unsigned int>(dof_indices_F.size());

    ke.resize (n_dofs , n_dofs);
    ke.zero();

    ke_Z.resize (n_dofs , n_dofs);
    ke_Z.zero();

    re.resize (n_dofs );
    re.zero();

  Real localPermeability=ComputeMaterialProprties(elem);
    
  permeability.assign( qrule->n_points() ,localPermeability );
    
  for (unsigned int i=0; i<phi.size(); i++)
    {
      for (unsigned int j=0; j<phi.size(); j++)
      {
        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        {
          ke(i,j) +=  JxW[qp] * ( dphi[j][qp] * ( permeability.at(qp) *  dphi[i][qp] ) );
        }
      }
    }

  for (auto side : elem->side_index_range())
  {
      if (elem->neighbor_ptr(side) == nullptr)
      {
        fe_face_F->reinit(elem, side);
        for(int k =0; k < _neumannSidesetIds.size(); k++)
        {
          if ( _meshBase[0].get_boundary_info().has_boundary_id (elem, side, _neumannSidesetIds.at(k)) )
          {
            for (unsigned int qp=0; qp<qface->n_points(); qp++)
              for (unsigned int i=0; i != n_dofs; i++){

               // std::cout<<" _neumannSidesetIds.at("<<k<<")"<<"is "<< _neumannSidesetIds.at(k)<<std::endl;
                
                re(i) += JxW_face[qp] * _neumannFunctions.at(k)[0].value(0.0,q_points_face[qp]) * phi_face[i][qp];
              }
          }
        }
      }
  }

  _dof_mapF[0].constrain_element_matrix_and_vector (ke, re, dof_indices_F,false);
  _dof_mapF[0].constrain_element_matrix (ke_Z, dof_indices_FZ,false);

  if (ke_Z.m()>n_dofs)
  {
    for (int ii=0; ii<ke_Z.m(); ++ii)
      if (ke_Z(ii,ii)>0.5)
      {
        ke(ii,ii)=0.0;
      }
  }

  for (auto side : elem->side_index_range())
  {
    Elem const * otherElem=elem->neighbor_ptr(side);
    if (otherElem != nullptr)
    {
      int otherID=otherElem->subdomain_id();
      if (elemId!=otherID)
      {
        //std::cout<<"close to the boundary\n";
        //fe_face->reinit(elem, side);

        std::unique_ptr<const Elem> test=elem->build_side_ptr (side);
        Elem const * ciaoF=test.get();

        _dof_mapF[0].dof_indices(ciaoF, dof_indices_FB);

        std::unique_ptr<FEBase> fel  (FEBase::build(_dim_F-1, fe_type));
        std::unique_ptr<QBase> qrulel( QBase::build (QGAUSS,_dim_F-1,THIRD));
        fel->attach_quadrature_rule (qrulel.get());

        std::vector<Real>                       const & JxWl = fel->get_JxW();
        std::vector<std::vector<Real> >         const & phil = fel->get_phi();
        std::vector<Point>                      const & q_pointsl = fel->get_xyz();

        fel->reinit(ciaoF);
        
        DenseMatrix<Number> kel;
        DenseVector<Number> rel;
        kel.resize(phil.size(),phil.size());
        rel.resize(phil.size());
        kel.zero();
        rel.zero();

        for (int i=0; i<phil.size(); ++i)
          for (int j=0; j<phil.size(); ++j)
            for (int qp=0; qp<phil.at(0).size(); ++qp)
              rel(i)+=JxWl.at(qp)*phil.at(i).at(qp)*phil.at(j).at(qp);

        _dof_mapF[0].constrain_element_matrix_and_vector (kel,rel, dof_indices_FB,true);

        for (int i=0; i<kel.m();++i)
        {
          Real & aa=kel(i,i);
          if (aa<0.5)
            aa=1.0;
          else
          {
            aa=0.0;
            for (int j=0; j<kel.n(); ++j)
            {
              Real & bb=kel(i,j);
              if (bb<-1e-10)
              {
                bb=-1.0*bb;
              }
            }
          }
        }

        if (elemId==_block_0_id)
        {
          iF0.add_vector (rel  , dof_indices_FB  );
        }
        else if(elemId==_block_1_id)
        {
          iF1.add_vector (rel  , dof_indices_FB  );
        }
        else if (elemId==_block_2_id)
        {
          iF2.add_vector (rel  , dof_indices_FB  );
        }
        else
        {
          mooseError("elemId does not belong to any block");
        }



      }
        
      }
  }


  if (elemId ==_block_0_id)
  {
    count_0=count_0+1;
    AF0.add_matrix (ke  , dof_indices_F  );
    bF0.add_vector (re  , dof_indices_F  );
  }
  else if (elemId ==_block_1_id)
  {
    count_1 = count_1+1;
    //std::cout<<"Here I am"<<std::endl;
    AF1.add_matrix (ke  , dof_indices_F  );
    bF1.add_vector (re  , dof_indices_F  );
  }
  else if (_has_block_2_id == true && elemId ==_block_2_id)
  {
    count_2 = count_2+1;
    AF2.add_matrix (ke  , dof_indices_F  );
    bF2.add_vector (re  , dof_indices_F  );
  }
  else
  {
    mooseError("elemId does not belong to any block");
  }



  // else
  // {
  //   mooseError("elemId does not belong to any block");
  // }

  }// elements loop

  AF0.close();
  AF1.close();
  AF2.close();
  bF0.close();
  bF1.close();
  bF2.close();
  iF0.close();
  iF1.close();
  iF2.close();

  for (int i=iF0.first_local_index(); i<iF0.last_local_index(); ++i)
  {
    if (iF0(i)>1e-10)
    {
      iF0.set( i,1./iF0(i) );
    }
    if (iF1(i)>1e-10)
    {
      iF1.set(i,1./iF1(i));
    }
    if (iF2(i)>1e-10)
    {
      iF2.set(i,1./iF2(i));
    }

  }

  iF0.close();
  iF1.close();
  iF2.close();

  std::cout << "ComputeFluxUtopia::assemble Fracture stop\n";
}



Real
ComputeFluxUtopia::ComputeMaterialProprties(const Elem *elem)
{
    Real permeability=0.0;

    for(int ll=0; ll<_vector_p.size(); ll++)
    {
      //std::cout<<elem->subdomain_id()<<std::endl;

        if (elem->subdomain_id()==_vector_p[ll])
        {
            permeability = _vector_value[ll];
        }
    }
    // std::cout<<permeability<<std::endl;
    return permeability;
}


