//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeFlux.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseMeshUtils.h"

// LibMesh includes

#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"

#include "helpers.h"

registerMooseObject("parrot2App", ComputeFlux);

defineLegacyParams(ComputeFlux);

InputParameters ComputeFlux::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<std::string>("material_name","material_name");

  params.addRequiredParam<int>("block_0_id", "block_0_id");

  params.addRequiredParam<int>("block_1_id", "block_1_id");
      
  params.addRequiredParam< std::vector<BoundaryName> >("dirichlet_nodeset_names","The name of the nodeset to create");

  params.addRequiredParam< std::vector<FunctionName> >("dirichlet_function_names", "The function.");

  params.addParam< std::vector<BoundaryName> >("neumann_sideset_names",std::vector<BoundaryName>(),"The ids of the nodeset to create");

  params.addParam< std::vector<FunctionName> >("neumann_function_names",std::vector<FunctionName>(), "The function.");

  params.addParam<std::string>("exodus_file", "the file name of the output");

  params.addParam<std::string>("nemesis_file", "the file name of the output");

  params.addParam<std::string>("variable_name_1", "the variable name of the output");

  params.addParam<std::string>("variable_name_2", "the variable name of the output");

  return params;

}

ComputeFlux::ComputeFlux(const InputParameters & parameters) :
GeneralUserObject(parameters),
_material_name ( getParam<std::string>("material_name") ),
_qrule(_assembly.qRule()),
_block_0_id( getParam< int >( "block_0_id" ) ),
_block_1_id( getParam< int >( "block_1_id" ) ),
_equationSystemsM( _fe_problem.es() ),
_dirichletNodesetNames ( getParam<std::vector<BoundaryName> >("dirichlet_nodeset_names" ) ),
_dirichletFunctionNames( getParam<std::vector<FunctionName> >("dirichlet_function_names") ),
_neumannSidesetNames   ( getParam<std::vector<BoundaryName> >("neumann_sideset_names"   ) ),
_neumannFunctionNames  ( getParam<std::vector<FunctionName> >("neumann_function_names"  ) ),
_has_exodus_file(  isParamValid("exodus_file")  ),
    // _has_nemesis_file( isParamValid("nemesis_file") ),
_has_variable_name_1( isParamValid("variable_name_1") ),
_has_variable_name_2( isParamValid("variable_name_2") )
{
  _sys_name           ="Flux_system";
  _var_name           ="Dummy_variable";
  _mat_domain_name_0  ="matrix_domain_0";
  _mat_domain_name_1  ="matrix_domain_1";
  _mat_binterp_name_0 ="bound_interp_0";
  _mat_binterp_name_1 ="bound_interp_1";
  _vec_flux_name_0    ="vector_flux_0";
  _vec_flux_name_1    ="vector_flux_1";
  _vec_rhs_name_0     ="vector_rhs_0";
  _vec_rhs_name_1     ="vector_rhs_1";
  _vec_boundary_name_0="vec_boundary_0";
  _vec_boundary_name_1="vec_boundary_1";
  _vec_dir            ="vector_dirichlet";

      if (_has_exodus_file)
        _exodus_filename=getParam<std::string>("exodus_file");
      
      // if (_has_nemesis_file)
      //   _nemesis_filename=getParam<std::string>("nemesis_file");
      
      if (_has_variable_name_1)
        _variable_name_1=getParam<std::string>("variable_name_1");

      if (_has_variable_name_2)
        _variable_name_2=getParam<std::string>("variable_name_2");

      if ( _dirichletNodesetNames.size() != _dirichletFunctionNames.size() )
      {
        mooseError("_dirichletNodesetNames.size() != _dirichletFunctionNames.size()");
      }

      if ( _neumannSidesetNames.size() != _neumannFunctionNames.size() )
      {
        mooseError("_neumannSidesetNames.size() != _neumannFunctionNames.size()");
      }
}

void ComputeFlux::init()
{
  std::cout<<"ComputeFlux::init() start\n";

  _meshBase=&_equationSystemsM.get_mesh();
  _linearImplicitSystemM = &_equationSystemsM.add_system<LinearImplicitSystem> (_sys_name.c_str());
  _v_var = _linearImplicitSystemM[0].add_variable (_var_name.c_str(), FIRST);
      
  _linearImplicitSystemM[0].add_matrix(_mat_domain_name_0.c_str());
  _linearImplicitSystemM[0].add_matrix(_mat_domain_name_1.c_str());
  _linearImplicitSystemM[0].add_matrix(_mat_binterp_name_0.c_str());
  _linearImplicitSystemM[0].add_matrix(_mat_binterp_name_1.c_str());

  // SparseMatrix<Number> & I=_linearImplicitSystemP[0].add_matrix("interpolator");

  _linearImplicitSystemM[0].add_vector(_vec_flux_name_0.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_flux_name_1.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_rhs_name_0.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_rhs_name_1.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_boundary_name_0.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_boundary_name_1.c_str());
  _linearImplicitSystemM[0].add_vector(_vec_dir.c_str());

  // // after this, the additional matrices are initialized
  std::cout<<"Before reinit, it may take a while...\n";
  _equationSystemsM.reinit();
  std::cout<<"Done\n";
  _equationSystemsM.print_info();

  _materialBase = &getMaterialByName(_material_name.c_str(),true);
  _flowAndTransport = dynamic_cast<FlowAndTransport *>(_materialBase);

  _dim = _meshBase[0].mesh_dimension();

  _dof_map = &_linearImplicitSystemM[0].get_dof_map();

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

  std::cout<<"ComputeFlux::init() stop\n";
}

void ComputeFlux::setDirichlet()
{
  std::cout<<"ComputeFlux::setDirichlet() start\n";

  NumericVector<Number> & d = _linearImplicitSystemM[0].get_vector(_vec_dir.c_str());
  d=1.0;

  // NumericVector<Number> & s=_linearImplicitSystemP[0].get_vector("solution_temp");
  // NumericVector<Number> & b=_linearImplicitSystemP[0].get_vector("rhs");
  // NumericVector<Number> & d=_linearImplicitSystemP[0].get_vector("dir");

  BoundaryInfo const & boundary_info = _meshBase[0].get_boundary_info();
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

        d.set(index,0.0);
        break;
      //       Point point=node[0];
      //       Real  value=_dirichletFunctions.at(i)[0].value(0.0,point) ;
            
      //       d.set(index,value);
      //       s.set(index,value);
      //       b.set(index,0.0);
      }
    }
  }

  d.close();
  // s.close();
  // b.close();

  std::cout<<"ComputeFlux::setDirichlet() stop\n";
}

void ComputeFlux::solve()
{
  std::cout<<"ComputeFlux::solve() start\n";

  SparseMatrix<Number> & A0 = _linearImplicitSystemM[0].get_matrix(_mat_domain_name_0.c_str());
  SparseMatrix<Number> & A1 = _linearImplicitSystemM[0].get_matrix(_mat_domain_name_1.c_str());

  SparseMatrix<Number> & I0 = _linearImplicitSystemM[0].get_matrix(_mat_binterp_name_0.c_str());
  SparseMatrix<Number> & I1 = _linearImplicitSystemM[0].get_matrix(_mat_binterp_name_1.c_str());

  NumericVector<Number> & f0 = _linearImplicitSystemM[0].get_vector(_vec_flux_name_0.c_str());
  NumericVector<Number> & f1 = _linearImplicitSystemM[0].get_vector(_vec_flux_name_1.c_str());
  NumericVector<Number> & b0 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_0.c_str());
  NumericVector<Number> & b1 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_1.c_str());
  NumericVector<Number> & d  = _linearImplicitSystemM[0].get_vector(_vec_dir.c_str());

  NumericVector<Number> & i0 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_0.c_str());
  NumericVector<Number> & i1 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_1.c_str());


  LinearImplicitSystem & linearImplicitSystem = _equationSystemsM.get_system<LinearImplicitSystem> ("nl0");

  std::unique_ptr<NumericVector<Number>> & solutionPointer(linearImplicitSystem.solution);
  NumericVector<Number> * sol=solutionPointer.get();

  A0.vector_mult(f0,sol[0]);
  f0*=-1.0;
  f0+=b0;
  f0*=d;
  f0.close();
  A1.vector_mult(f1,sol[0]);
  f1*=-1.0;
  f1+=b1;
  f1*=d;
  f1.close();

  std::cout << std::fixed;
  std::cout<<std::setprecision(6)<<f0.sum()<<std::endl;
  std::cout<<std::setprecision(6)<<f1.sum()<<std::endl;

  std::cout<<i0.sum()<<std::endl;
  std::cout<<i1.sum()<<std::endl;

  d=f0;
  d*=i0;
  I0.vector_mult(b0,d);

  b0.close();
  // sol[0]=b0;
  // sol[0].close();

  d=f1;
  d*=i1;
  I1.vector_mult(b1,d);

  b1.close();

  std::cout<<"ComputeFlux::solve() stop\n";
}

    void ComputeFlux::write()
    {
      std::cout<<"ComputeFlux::write() start\n";

      if(_has_exodus_file)
      {
        ExodusII_IO (_equationSystemsM.get_mesh()).write_equation_systems(_exodus_filename.c_str(), _equationSystemsM); //output_filename.c_str()
      }

      // if(_has_nemesis_file)
      // {
      //   Nemesis_IO (_equationSystemsP[0].get_mesh()).write_equation_systems(_nemesis_filename.c_str(), _equationSystemsP[0]); //output_filename.c_str()
      // }

      std::cout<<"ComputeFlux::write() stop\n";
    }

    void ComputeFlux::del()
    {
      std::cout<<"ComputeFlux::del() start\n";

      // //_equationSystemsP[0].delete_system(_sys_name.c_str());
      // _equationSystemsP[0].clear();

      std::cout<<"ComputeFlux::del() stop\n";  
    }

    void ComputeFlux::initialize()
    {

      init();
      assemble();
      setDirichlet();
      solve();
      write();
      
      //_equationSystemsM.reinit();
      
      if(_has_variable_name_1 && _has_variable_name_2)

      {
        //copyVariableFl(_equationSystemsM, _var_name, _equationSystemsM, _variable_name);

        NumericVector<Number> & b0 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_0.c_str());
        NumericVector<Number> & b1 = _linearImplicitSystemM[0].get_vector(_vec_rhs_name_1.c_str());

        //NumericVector<Number> & f0 = _linearImplicitSystemM[0].get_vector(_vec_flux_name_0.c_str());
        //NumericVector<Number> & f1 = _linearImplicitSystemM[0].get_vector(_vec_flux_name_1.c_str());

        MooseVariableFEBase  & _flux_var_1 = _fe_problem.getVariable(0, _variable_name_1, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
        MooseVariableFEBase  & _flux_var_2 = _fe_problem.getVariable(0, _variable_name_2, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
        MooseVariableFEBase  & sol_var = _fe_problem.getVariable(0, "pressure", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

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

void ComputeFlux::assemble()
{
  std::cout << "ComputeFlux::assemble start\n";

  SparseMatrix <Number> & A0=_linearImplicitSystemM[0].get_matrix( _mat_domain_name_0.c_str() );
  SparseMatrix <Number> & A1=_linearImplicitSystemM[0].get_matrix( _mat_domain_name_1.c_str() );
  SparseMatrix<Number> & I0 = _linearImplicitSystemM[0].get_matrix(_mat_binterp_name_0.c_str());
  SparseMatrix<Number> & I1 = _linearImplicitSystemM[0].get_matrix(_mat_binterp_name_1.c_str());

  NumericVector<Number> & b0=_linearImplicitSystemM[0].get_vector( _vec_rhs_name_0.c_str() );
  NumericVector<Number> & b1=_linearImplicitSystemM[0].get_vector( _vec_rhs_name_1.c_str() );

  NumericVector<Number> & i0 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_0.c_str());
  NumericVector<Number> & i1 = _linearImplicitSystemM[0].get_vector(_vec_boundary_name_1.c_str());

  A0.zero();
  A1.zero();
  b0.zero();
  b1.zero();
  i0.zero();
  i1.zero();

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

  for ( ; el != end_el; ++el)
  {
    Elem const * elem = *el;
    fe->reinit (elem);

    int elemId=elem->subdomain_id();

    _dof_map[0].dof_indices(elem, dof_indices);
    _dof_map[0].dof_indices(elem, dof_indices_Z);
  //   dof_indices_I=dof_indices;

    unsigned int const n_dofs = cast_int<unsigned int>(dof_indices.size());

    ke.resize (n_dofs , n_dofs);
    ke.zero();

    ke_Z.resize (n_dofs , n_dofs);
    ke_Z.zero();

    re.resize (n_dofs );
    re.zero();
  //   ke_I.resize (n_dofs , n_dofs);
  //   ke_I.zero();

    _flowAndTransport[0].getPermeability( q_points,permeability);
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
              for (unsigned int i=0; i != n_dofs; i++)
                re(i) += JxW_face[qp] * _neumannFunctions.at(k)[0].value(0.0,q_points_face[qp]) * phi_face[i][qp];
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

        std::unique_ptr<FEBase> fel  (FEBase::build(_dim-1, fe_type));
        std::unique_ptr<QBase> qrulel( QBase::build (_qrule->type(),_dim-1,_qrule->get_order()));
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
          i0.add_vector (rel  , dof_indices_B  );
          for (int i=0; i<kel.m();++i)
            for (int j=0; j<kel.n();++j)
              I0.set(dof_indices_B.at(i),dof_indices_B.at(j),kel(i,j));
        }
        else if (elemId==_block_1_id)
        {
          i1.add_vector (rel  , dof_indices_B  );
          for (int i=0; i<kel.m();++i)
            for (int j=0; j<kel.n();++j)
              I1.set(dof_indices_B.at(i),dof_indices_B.at(j),kel(i,j));
        }
        else
        {
          mooseError("elemId does not belong to any block");
        }

      }
        
      }
  }

  if (elemId==_block_0_id)
  {
    A0.add_matrix (ke  , dof_indices  );
    b0.add_vector (re  , dof_indices  );
  }
  else if (elemId==_block_1_id)
  {
    A1.add_matrix (ke  , dof_indices  );
    b1.add_vector (re  , dof_indices  );
  }
  else
  {
    mooseError("elemId does not belong to any block");
  }

  }// elements loop

  A0.close();
  A1.close();
  I0.close();
  I1.close();
  b0.close();
  b1.close();
  i0.close();
  i1.close();

  for (int i=i0.first_local_index(); i<i0.last_local_index(); ++i)
  {
    if (i0(i)>1e-10)
    {
      i0.set( i,1./i0(i) );
      //std::cout<<1./i0(i)<<std::endl;
    }
    if (i1(i)>1e-10)
    {
      i1.set(i,1./i1(i));
    }
  }

  i0.close();
  i1.close();

  std::cout << "ComputeFlux::assemble stop\n";
}


