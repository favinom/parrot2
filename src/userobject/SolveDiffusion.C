//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolveDiffusion.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseMeshUtils.h"

// LibMesh includes

#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"


#include "helpers.h"



// inline void convert_vec_SD(libMesh::NumericVector<libMesh::Number> &lm_vec,  utopia::UVector &utopia_vec) {
//     using namespace libMesh;
//     Vec p_vec = cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&lm_vec)->vec();
//     utopia_vec.wrap(p_vec);
// }

// inline void convert_mat_SD(libMesh::SparseMatrix<libMesh::Number> &lm_mat, utopia::USparseMatrix &utopia_mat) {
//     using namespace libMesh;

//     Mat p_mat = cast_ptr<libMesh::PetscMatrix<libMesh::Number> *>(&lm_mat)->mat();
//     utopia::convert(p_mat, utopia_mat);
//     }


registerMooseObject("parrot2App", SolveDiffusion);

defineLegacyParams(SolveDiffusion);

InputParameters SolveDiffusion::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<std::string>("material_name","material_name");
  
  params.addRequiredParam< std::vector<BoundaryName> >("dirichlet_nodeset_names","The name of the nodeset to create");

  params.addRequiredParam< std::vector<FunctionName> >("dirichlet_function_names", "The function.");

  params.addParam< std::vector<BoundaryName> >("neumann_sideset_names",std::vector<BoundaryName>(),"The ids of the nodeset to create");

  params.addParam< std::vector<FunctionName> >("neumann_function_names",std::vector<FunctionName>(), "The function.");

  params.addRequiredParam<int>("solver_type","solver_type");

  params.addParam<std::string>("exodus_file", "the file name of the output");

  params.addParam<std::string>("nemesis_file", "the file name of the output");

  params.addParam<std::string>("variable_name", "the variable name of the output");

  params.addParam<bool>("multisubdomain",false,"multisubdomain");

  return params;
}

SolveDiffusion::SolveDiffusion(const InputParameters & parameters) :
GeneralUserObject(parameters),
_material_name ( getParam<std::string>("material_name") ),
_qrule(_assembly.qRule()),
_dirichletNodesetNames ( getParam<std::vector<BoundaryName> >("dirichlet_nodeset_names" ) ),
_dirichletFunctionNames( getParam<std::vector<FunctionName> >("dirichlet_function_names") ),
_neumannSidesetNames   ( getParam<std::vector<BoundaryName> >("neumann_sideset_names"   ) ),
_neumannFunctionNames  ( getParam<std::vector<FunctionName> >("neumann_function_names"  ) ),
_solverType(getParam<int>("solver_type")),
_has_exodus_file(  isParamValid("exodus_file")  ),
_has_nemesis_file( isParamValid("nemesis_file") ),
_has_variable_name( isParamValid("variable_name") ),
_multisubdomain(getParam<bool>("multisubdomain"))
/*
_aux_var_names(getParam<std::vector<AuxVariableName>>("aux_variable")),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
_boundary_D_ids(getParam<std::vector<boundary_id_type>>("boundary_D_bc")),
_boundary_N_ids(getParam<std::vector<boundary_id_type>>("boundary_N_bc")),
_value_N_bc(getParam<std::vector<Real>>("value_N_bc")),
_value_D_bc(getParam<std::vector<Real>>("value_D_bc")),

_hasMeshModifier( isParamValid("fractureMeshModifier") ),
_conservativeScheme( getParam<bool>("conservative") ),
*/
{
  _sys_name="Diffusion";
  _var_name="pressure";

  if (_has_exodus_file)
    _exodus_filename=getParam<std::string>("exodus_file");
  if (_has_nemesis_file)
    _nemesis_filename=getParam<std::string>("nemesis_file");
  if (_has_variable_name)
    _variable_name=getParam<std::string>("variable_name");

  if ( _dirichletNodesetNames.size() != _dirichletFunctionNames.size() )
  {
    mooseError("_dirichletNodesetNames.size() != _dirichletFunctionNames.size()");
  }

  if ( _neumannSidesetNames.size() != _neumannFunctionNames.size() )
  {
    mooseError("_neumannSidesetNames.size() != _neumannFunctionNames.size()");
  }
}

  // Parallel::Communicator const & pp_comm=meshBaseM.comm();
  // _mesh=new DistributedMesh(pp_comm);
  // _mesh[0].read("mesh.xda");
  // _meshBase=_mesh;
//  MeshBase & meshBaseM=_equationSystemsM[0].get_mesh();


void SolveDiffusion::init()
{
  std::cout<<"SolveDiffusion::init() start\n";

  _equationSystemsM=&_fe_problem.es();
  _meshBase=&_equationSystemsM[0].get_mesh();
  
  _pp_comm.push_back( &_meshBase[0].comm() );

  _equationSystemsP=new EquationSystems(_meshBase[0]);
  _linearImplicitSystemP = &_equationSystemsP[0].add_system<LinearImplicitSystem> (_sys_name.c_str());
  _p_var = _linearImplicitSystemP[0].add_variable (_var_name.c_str(), FIRST);
  
  SparseMatrix<Number> & I=_linearImplicitSystemP[0].add_matrix("interpolator");

  _linearImplicitSystemP[0].add_vector("solution_temp");
  _linearImplicitSystemP[0].add_vector("rhs");
  _linearImplicitSystemP[0].add_vector("dir");

  // after this, the additional matrices are initialized
  _equationSystemsP[0].init();
  _equationSystemsP[0].print_info();

  _materialBase = &getMaterialByName(_material_name.c_str(),true);
  _flowAndTransport = dynamic_cast<FlowAndTransport *>(_materialBase);

  _dim = _meshBase[0].mesh_dimension();

  _dof_map = &_linearImplicitSystemP[0].get_dof_map();

  _parrotSolver=new ParrotSolver( _solverType );

  for (int i=0; i<_dirichletFunctionNames.size(); ++i)
  {
    Function const & func=getFunctionByName( _dirichletFunctionNames.at(i) );
    _dirichletFunctions.push_back(&func);
  }

  for (int i=0; i<_neumannFunctionNames.size(); ++i)
  {
    Function const & func=getFunctionByName( _neumannFunctionNames.at(i) );
    _neumannFunctions.push_back(&func);
  }

  _dirichletSidesetIds = MooseMeshUtils::getBoundaryIDs(_meshBase[0], _dirichletNodesetNames, true);
  _neumannSidesetIds   = MooseMeshUtils::getBoundaryIDs(_meshBase[0], _neumannSidesetNames,   true);

  std::cout<<"SolveDiffusion::init() stop\n";
}


void SolveDiffusion::setDirichlet()
{
  std::cout<<"SolveDiffusion::setDirichlet() start\n";

  NumericVector<Number> & s=_linearImplicitSystemP[0].get_vector("solution_temp");
  NumericVector<Number> & b=_linearImplicitSystemP[0].get_vector("rhs");
  NumericVector<Number> & d=_linearImplicitSystemP[0].get_vector("dir");

  BoundaryInfo & boundary_info = _meshBase[0].get_boundary_info();

  for ( auto & node : as_range(_meshBase->active_nodes_begin(), _meshBase->active_nodes_end()) )
  {
    for (int i=0; i<_dirichletSidesetIds.size(); ++i)
    {
      if ( boundary_info.has_boundary_id(node,_dirichletSidesetIds.at(i)) )
      {
        std::vector<dof_id_type> di;
        _dof_map[0].dof_indices (node,di, 0);
        if (di.size()!=1)
        {
          mooseError("di.size()!=1");
        }
        dof_id_type index=di.at(0);
        _dirIds.push_back(index);

        Point point=node[0];
        Real  value=_dirichletFunctions.at(i)[0].value(0.0,point) ;
        
        d.set(index,value);
        s.set(index,value);
        b.set(index,0.0);
        break;
      }
    }
  }

  d.close();
  s.close();
  b.close();

  std::cout<<"SolveDiffusion::setDirichlet() stop\n";

}



void SolveDiffusion::solve()
{
  std::cout<<"SolveDiffusion::solve() start\n";

  SparseMatrix<Number> & M=_linearImplicitSystemP[0].get_system_matrix();
  SparseMatrix<Number> & I=_linearImplicitSystemP[0].get_matrix("interpolator");

  NumericVector<Number> & s=_linearImplicitSystemP[0].get_vector("solution_temp");
  NumericVector<Number> & b=_linearImplicitSystemP[0].get_vector("rhs");
  NumericVector<Number> & d=_linearImplicitSystemP[0].get_vector("dir");
  b.add(d);
  
  M.zero_rows(_dirIds,1.0);

  SparseMatrix<Number>  * M_ptr = &M;
  NumericVector<Number> * b_ptr = &b;
  NumericVector<Number> * s_ptr = &s;

  std::unique_ptr<NumericVector<Number>> & solutionPointer(_linearImplicitSystemP[0].solution);
  NumericVector<Number> * sol=solutionPointer.get();

  _parrotSolver->setMatrixAndVectors(M_ptr,b_ptr,s_ptr);
  _parrotSolver->setConstantMatrix(true);
  _parrotSolver->solve();

  I.vector_mult(sol[0],s);
  
  sol[0].close();
  _linearImplicitSystemP[0].update();

  // utopia::UVector U_b;
  // utopia::USparseMatrix U_m;

  // convert_vec_SD(b,U_b);
  // convert_mat_SD(M,U_m);


  // U_m.write("system_matrix.m");
  // U_b.write("system_vector.m");

  std::cout<<"SolveDiffusion::solve() stop\n";
}

void SolveDiffusion::write()
{
  std::cout<<"SolveDiffusion::write() start\n";

  if(_has_exodus_file)
  {
    ExodusII_IO (_equationSystemsP[0].get_mesh()).write_equation_systems(_exodus_filename.c_str(), _equationSystemsP[0]); //output_filename.c_str()
  }

  if(_has_nemesis_file)
  {
    Nemesis_IO (_equationSystemsP[0].get_mesh()).write_equation_systems(_nemesis_filename.c_str(), _equationSystemsP[0]); //output_filename.c_str()
  }

  std::cout<<"SolveDiffusion::write() stop\n";
}

void SolveDiffusion::del()
{
  std::cout<<"SolveDiffusion::del() start\n";

  //_equationSystemsP[0].delete_system(_sys_name.c_str());
  _equationSystemsP[0].clear();


  delete _parrotSolver;
  delete _equationSystemsP;

  std::cout<<"SolveDiffusion::del() stop\n";  
}

void SolveDiffusion::initialize()
{

  init();
  assemble();
  setDirichlet();
  solve();
  write();
  _equationSystemsM[0].reinit();
  if(_has_variable_name)
  {
    copyVariable(_equationSystemsP[0], _var_name, _equationSystemsM[0], _variable_name);
  }
  del();

  std::cout<<"END initialize\n";
}

void SolveDiffusion::assemble()
{
  std::cout << "SolveDiffusion::assemble start\n";

  SparseMatrix<Number> & M=_linearImplicitSystemP[0].get_system_matrix();
  SparseMatrix<Number> & I=_linearImplicitSystemP[0].get_matrix("interpolator");

  NumericVector<Number> & b=_linearImplicitSystemP[0].get_vector("rhs");

  M.zero();
  I.zero();

  FEType fe_type = _dof_map[0].variable_type(0);
  
  std::unique_ptr<FEBase> fe (FEBase::build(_dim, fe_type));
  std::unique_ptr<QBase> qrule( QBase::build (_qrule->type(),_dim,_qrule->get_order()));
  fe->attach_quadrature_rule (qrule.get());

  std::unique_ptr<FEBase> fe_face (FEBase::build(_dim, fe_type));
  std::unique_ptr<QBase> qface(fe_type.default_quadrature_rule(_dim-1));
  fe_face->attach_quadrature_rule (qface.get());

  std::vector<Real>                       const & JxW = fe->get_JxW();
  std::vector<std::vector<Real> >         const & phi = fe->get_phi();
  std::vector<std::vector<RealGradient> > const & dphi = fe->get_dphi();
  std::vector<Point>                      const & q_points = fe->get_xyz();

  std::vector<Real>                       const & JxW_face      = fe_face->get_JxW();
  std::vector<std::vector<Real> >         const & phi_face      = fe_face->get_phi();
  std::vector<std::vector<RealGradient> > const & dphi_face     = fe_face->get_dphi();
  std::vector<Point>                      const & q_points_face = fe_face->get_xyz();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_I;

  DenseMatrix<Number> ke;
  DenseMatrix<Number> ke_I;
  DenseVector<Number> re;

  MeshBase::const_element_iterator           el = _meshBase[0].active_local_elements_begin();
  MeshBase::const_element_iterator const end_el = _meshBase[0].active_local_elements_end();

  std::vector<Number> permeability;

  for ( ; el != end_el; ++el)
  {
    Elem const * elem = *el;
    
    fe->reinit (elem);
    _dof_map[0].dof_indices(elem, dof_indices);
    dof_indices_I=dof_indices;

    unsigned int const n_dofs = cast_int<unsigned int>(dof_indices.size());

    ke.resize (n_dofs , n_dofs);
    ke.zero();
    re.resize (n_dofs );
    re.zero();
    ke_I.resize (n_dofs , n_dofs);
    ke_I.zero();


    int subdomain_id = elem->subdomain_id();

    
    if (_multisubdomain){
      _flowAndTransport[0].getPermeabilityId(q_points,permeability,subdomain_id);
    }
    else{
      _flowAndTransport[0].getPermeability(q_points,permeability);
    }

    for (unsigned int i=0; i<phi.size(); i++)
      for (unsigned int j=0; j<phi.size(); j++)
        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        {
          ke(i,j) +=  JxW[qp] * ( dphi[j][qp] * ( permeability.at(qp) *  dphi[i][qp] ) );
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
              for (unsigned int i=0; i != n_dofs; i++)
                re(i) += JxW_face[qp] * _neumannFunctions.at(k)[0].value(0.0,q_points_face[qp]) * phi_face[i][qp];
          }
        }
      }
    }

    _dof_map[0].constrain_element_matrix_and_vector (ke, re, dof_indices,false);
    _dof_map[0].constrain_element_matrix (ke_I, dof_indices_I);

    for (int i=0; i<ke_I.m();++i)
    {
      Real & aa=ke_I(i,i);
      if (aa<0.5)
        aa=1.0;
      else
      {
        aa=0.0;
        for (int j=0; j<ke_I.n(); ++j)
        {
          Real & bb=ke_I(i,j);
          if (bb<-1e-10)
          {
            bb=-1.0*bb;
          }
        }
      }
    }

    M.add_matrix (ke  , dof_indices  );
    //I.add_matrix (ke_I, dof_indices_I);
    for (int i=0; i<ke_I.m();++i)
      for (int j=0; j<ke_I.n();++j)
        I.set(dof_indices_I.at(i),dof_indices_I.at(j),ke_I(i,j));

    b.add_vector (re, dof_indices);

  }// elements loop

  M.close();
  I.close();
  b.close();




  std::cout << "SolveDiffusion::assemble stop\n";
}
