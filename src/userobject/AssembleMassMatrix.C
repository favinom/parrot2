//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssembleMassMatrix.h"

#include "FEProblem.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "Assembly.h"
// #include "MooseVariableFEBase.h"

#include "libmesh/quadrature_gauss.h"

registerMooseObject("parrot2App", AssembleMassMatrix);


template <>
InputParameters
validParams<AssembleMassMatrix>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addRequiredParam<std::string>("material_name","material_name");
  params.addRequiredParam<bool>("constrain_matrix","constrain_matrix");
  params.addParam<std::string>("dc_boundaries", "-1", "Dirichlet Boundary ID");
  params.addParam<std::string>("dc_variables" , "-1", "Variable to which given BC_id applies");
  params.addRequiredParam<std::vector<Real>>("value_D_bc", "The value of Dirichlet");
  params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators"); 
  return params;
}


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

AssembleMassMatrix::AssembleMassMatrix(const InputParameters & parameters) :
GeneralUserObject(parameters),
_material_name ( getParam<std::string>("material_name") ),
_userObjectName(getParam<UserObjectName>("operator_userobject")),
_constrainMatrices(getParam<bool>("constrain_matrix")),
_code_dof_map(true),
_qrule(_assembly.qRule()),
_dc_var(getParam<std::string>("dc_variables")),
_value_D_bc(getParam<std::vector<Real>>("value_D_bc"))


{

  std::vector<std::string> tmp = split_string(parameters.get<std::string>("dc_boundaries"), ' ');
    for(auto str_tmp=tmp.begin(); str_tmp != tmp.end(); str_tmp++)
    {
        _dc_boundary_id.push_back(atoi(str_tmp->c_str()));
    }
}

void AssembleMassMatrix::execute()
{
  assemble_mass_matrix();
  determine_dc_bnd_var_id(AssembleMassMatrix::split_string(_dc_var, ' '));

  find_boundary(_dc_boundary_id);
};

void AssembleMassMatrix::initialize()
{

  TransientNonlinearImplicitSystem const & _system2 = _fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
       
      // _interpolator->init();
    _equationSystemsT=&_fe_problem.es();

    _meshBase=&_equationSystemsT[0].get_mesh();

    _pp_comm.push_back( &_meshBase[0].comm() );

    _linearImplicitSystemT = &_equationSystemsT->add_system<LinearImplicitSystem> ("Transport");

    _test_var = _linearImplicitSystemT[0].add_variable ("test_var", FIRST);

    SparseMatrix<Number> & _interpolator = _linearImplicitSystemT[0].add_matrix("Interpolator");
    // SparseMatrix<Number> & _mass_matrix  = _linearImplicitSystemT[0].add_matrix("MassMatrix");
    // SparseMatrix<Number> & _poro_mass_matrix  = _linearImplicitSystemT[0].add_matrix("Poro_mass_matrix");

    // SparseMatrix<Number> & _lump_mass_matrix  = _linearImplicitSystemT[0].add_matrix("Lump_mass_matrix");
    SparseMatrix<Number> & _poro_lump_mass_matrix  = _linearImplicitSystemT[0].add_matrix("Poro_lump_mass_matrix");
    SparseMatrix<Number> & _hanging_interpolator  = _linearImplicitSystemT[0].add_matrix("Hanging Interpolator");
    SparseMatrix<Number> & _stab_matrix  =  _linearImplicitSystemT[0].add_matrix("_stab_matrix");

    _equationSystemsT[0].reinit();
    _equationSystemsT[0].print_info();

}

void AssembleMassMatrix::assemble_mass_matrix(){

    _console << "Assemble_Mass_matrix() begin "  << std::endl;

    _materialBase = &getMaterialByName(_material_name.c_str(),true);
    _flowAndTransport = dynamic_cast<FlowAndTransport *>(_materialBase);

    DofMap   const & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();

    StoreOperators & storeOperatorsUO=(_fe_problem.getUserObject<StoreOperators>(_userObjectName));
    _hanging_vec           = storeOperatorsUO.HangVec();
    _hanging_vec->init(dof_map.n_dofs(), dof_map.n_local_dofs());
    _hanging_vec->zero();
    _hanging_vec->add(1.0);

    _equationSystemsT=&_fe_problem.es();



    _linearImplicitSystemT = &_equationSystemsT->get_system<LinearImplicitSystem> ("Transport");



    SparseMatrix<Number> & _interpolator = _linearImplicitSystemT[0].get_matrix("Interpolator");
    // SparseMatrix<Number> & _mass_matrix  = _linearImplicitSystemT[0].get_matrix("MassMatrix");
    // SparseMatrix<Number> & _poro_mass_matrix  = _linearImplicitSystemT[0].get_matrix("Poro_mass_matrix");

    // SparseMatrix<Number> & _lump_mass_matrix  = _linearImplicitSystemT[0].get_matrix("Lump_mass_matrix");
    SparseMatrix<Number> & _poro_lump_mass_matrix  = _linearImplicitSystemT[0].get_matrix("Poro_lump_mass_matrix");
    SparseMatrix<Number> & _hanging_interpolator  = _linearImplicitSystemT[0].get_matrix("Hanging Interpolator");


    // Get a constant reference to the mesh object.
    MeshBase     const & mesh = _fe_problem.es().get_mesh();
    unsigned int const   dim  = mesh.mesh_dimension();

    // Get a reference to our system.
    TransientNonlinearImplicitSystem const & _system = _fe_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType const & fe_type = _system.get_dof_map().variable_type(0);
    //FEType fe_type = system.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    //QGauss qrule (dim, fe_type.default_quadrature_order());
    std::unique_ptr<QBase> qrule( QBase::build (_qrule->type(),dim,_qrule->get_order()));
    //QGauss qrule (dim, TENTH);
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (qrule.get());

    const std::vector<Real>& JxW      = fe->get_JxW();    
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    //const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
    const std::vector<Point>& q_points = fe->get_xyz();
    //const DofMap& dof_map = system.get_dof_map();

    // std::vector<dof_id_type> dof_indices;
    // std::vector<dof_id_type> dof_indices_l;
    // std::vector<dof_id_type> dof_indices_p;
    std::vector<dof_id_type> dof_indices_l_p;
    std::vector<dof_id_type> dof_indices_i;
    std::vector<dof_id_type> dof_indices_h;

    // DenseMatrix<Number> Me;
    // DenseMatrix<Number> Me_p;
    // DenseMatrix<Number> Me_l;
    DenseMatrix<Number> Me_l_p;
    DenseMatrix<Number> Me_i, Me_h;
  

    _interpolator.zero();
    // _mass_matrix.zero();
    // _poro_mass_matrix.zero();
    // _lump_mass_matrix.zero();
    _poro_lump_mass_matrix.zero();
    // _poro_mass_matrix.zero();

      


    MeshBase::const_element_iterator       el    = mesh.active_local_elements_begin();
    MeshBase::const_element_iterator const end_el = mesh.active_local_elements_end();


    std::vector<Real> poroVec;

    for ( ; el != end_el; ++el)
    {
    	const Elem * elem = *el;
    	fe->reinit (elem);

    	// dof_map.dof_indices(elem, dof_indices);
    	// dof_map.dof_indices(elem, dof_indices_l);
    	// dof_map.dof_indices(elem, dof_indices_p);
    	dof_map.dof_indices(elem, dof_indices_l_p);
    	dof_map.dof_indices(elem, dof_indices_i);
      dof_map.dof_indices(elem, dof_indices_h);

    	int const loc_n=dof_indices_l_p.size();


    	// Me.resize(loc_n,loc_n);
    	// Me_p.resize(loc_n,loc_n);
    	// Me_l.resize(loc_n,loc_n);
    	Me_l_p.resize(loc_n,loc_n);
    	Me_i.resize(loc_n,loc_n);
      Me_h.resize(loc_n,loc_n);

    	// Me.zero();
    	// Me_p.zero();
    	// Me_l.zero();
    	Me_l_p.zero();
    	Me_i.zero();
      Me_h.zero();

      std::vector<Number>  poroVec;
        
      _flowAndTransport[0].getPorosity(q_points,poroVec);
        

    	for (unsigned int i=0; i<phi.size(); i++)
    	{
    		for (unsigned int j=0; j<phi.size(); j++)
    		{
    			for (unsigned int qp=0; qp<qrule->n_points(); qp++)
    			{
    				// Me  (i,j)   +=  JxW[qp] * phi[i][qp] * phi[j][qp];
    				// Me_p(i,j)   += poroVec.at(qp) * JxW[qp] * phi[i][qp] * phi[j][qp];
        //     //std::cout<<"poroVec.at(qp)"<<poroVec.at(qp)<<std::endl;
    				// Me_l(i,i)   += JxW[qp] * phi[i][qp] * phi[j][qp];
    				Me_l_p(i,i) += poroVec.at(qp) * JxW[qp] * phi[i][qp] * phi[j][qp];
            //Me_h(i,i)    = 0.5;
    			}
    		}
    	}

    	if (_constrainMatrices)
    	{
    		// dof_map.constrain_element_matrix(Me,dof_indices);
    		// dof_map.constrain_element_matrix(Me_p,dof_indices_p);
    		// dof_map.constrain_element_matrix(Me_l,dof_indices_l);
    		dof_map.constrain_element_matrix(Me_l_p,dof_indices_l_p);
    	}
    	{
    		dof_map.constrain_element_matrix(Me_i,dof_indices_i);
        dof_map.constrain_element_matrix(Me_h, dof_indices_h);
      

        //std::cout<<"a"<<Me_h.m()<<"b"<<Me_h.n()<<"c"<<dof_indices_h.size()<<std::endl;

        for (int i=0; i<Me_h.m(); ++i)
        {

          for (int j=0; j<Me_h.n(); ++j)
           {
            if(Me_h(i,j)==1){

              _hanging_vec->set(dof_indices_h.at(i),0.0);

            }
            if (Me_h(i,j)<0){

              // std::cout<<"begin"<<std::endl;
              // std::cout<<"dof_indices.at(i)"<<dof_indices.at(i)<<std::endl;
              // std::cout<<"dof_indices.at(j)"<<dof_indices.at(j)<<std::endl;
              // std::cout<<"value"<<Me_h(i,j)<<std::endl;
              //std::cout<<"end"<<std::endl;
              //_hanging_vec->set(dof_indices.at(j),0.0);

              Real value = -1.0 * Me_h(i,j);
              
              _hanging_interpolator.set(dof_indices_h.at(i),dof_indices_h.at(j),value);

            }
          }
        }
      }

  
     {
    		
        for (int i=0; i<Me_i.m(); ++i)
    		{
          //Me_h(i,i)=1.0;

    			if (Me_i(i,i)<0.5)
    			{
    				Me_i(i,i)=1;
    			}
    			else
    			{
    				for (int j=0; j<Me_i.n(); ++j)
    				{
    					Me_i(i,j)=-1.0*Me_i(i,j);
    				}
    				Me_i(i,i)=0.0;
    			}
    		}
    	}

      

    	// (_mass_matrix).add_matrix(Me, dof_indices);
    	// (_poro_mass_matrix).add_matrix (Me_p, dof_indices_p);
    	// (_lump_mass_matrix).add_matrix (Me_l, dof_indices_l);
    	(_poro_lump_mass_matrix).add_matrix (Me_l_p, dof_indices_l_p);
    	(_interpolator).add_matrix (Me_i, dof_indices_i);

      //(*_hanging_interpolator).add_matrix (Me_h, dof_indices_h);
    }



   // (_mass_matrix).close();
   // (_poro_mass_matrix).close();
   (_poro_lump_mass_matrix).close();
   // (_lump_mass_matrix).close();
   (_interpolator).close();
   (_hanging_interpolator).close();

    _hanging_vec->close();

   
    _console << "Assemble_Mass_matrix() end "  << std::endl;


}




void
AssembleMassMatrix::find_boundary(std::vector<int> &_dc_boundary_id){




  _console << "AssembleMassMatrix::find_boundary begin "  << std::endl;


  ConstBndNodeRange & bnd_nodes = *_fe_problem.mesh().getBoundaryNodeRange();

  NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();

  const DofMap & dof_map = _nl.dofMap(); 

  unsigned int i = 0;

  StoreOperators & storeOperatorsUO=(_fe_problem.getUserObject<StoreOperators>(_userObjectName));

  _bc_vec                = storeOperatorsUO.BcVec();

  _bc_vec->init(dof_map.n_dofs(), dof_map.n_local_dofs());
  
  _bc_vec->zero();

  _bc_vec->add(1.0);

  _value_bc_vec                = storeOperatorsUO.ValueBcVec();

  _value_bc_vec->init(dof_map.n_dofs(), dof_map.n_local_dofs());
  
  _value_bc_vec->zero();



    // std::cout<<"_dc_variables_id"<< _dc_variables_id[0].size()<<std::endl;
    // std::cout<<"_dc_boundary_id"<< _dc_boundary_id.size()<<std::endl;
    for(auto boundary = _dc_boundary_id.begin(); boundary != _dc_boundary_id.end(); ++boundary, i++)
      {
        // iterate just over boundary nodes
            for (const auto & bnode : bnd_nodes)
            {
                  libMesh::Node * current_node = bnode->_node;

                  // check if node is in active boundary list
                  if (_fe_problem.mesh().isBoundaryNode(current_node->id(), *boundary))
                  {
                    // loop over all variables at this node

                    for (auto v = 0; v < _fe_problem.getNonlinearSystemBase().nVariables(); v++)
                    {
                      const Variable & var = _nl.system().variable(v);
                      
                      unsigned int var_num = var.number();
                        //std::cout<<"nnnnnnn"<< var_num <<std::endl;

                      // see if this variable has any dofs at this node
                      if (current_node->n_dofs(_fe_problem.getNonlinearSystemBase().number(), var_num) > 0)
                      {
                        // check if given variable has BC on node

                        if(std::find(_dc_variables_id[i].begin(), _dc_variables_id[i].end(), var_num) != _dc_variables_id[i].end())
                        {

                          _bc_vec->set(current_node->dof_number(_fe_problem.getNonlinearSystemBase().number(), var_num, 0), 0.0);
                          _value_bc_vec->set(current_node->dof_number(_fe_problem.getNonlinearSystemBase().number(), var_num, 0), _value_D_bc.at(0));
                        }
                    }
                }
            } 
        }
    }

    _bc_vec->close();


    _console << "AssembleMassMatrix::find_boundary end "  << std::endl;
  

     
     //auto it = std::find(zero_rows.begin(), zero_rows.end(), row);

    //std::cout<<"zero_rows"<< zero_rows.size()<<std::endl;
}


void 
AssembleMassMatrix::determine_dc_bnd_var_id(const std::vector<std::string> & BC_var){
    // automatic fill-in
     NonlinearSystemBase & _nl = _fe_problem.getNonlinearSystemBase();

    std::vector<int> vec(_nl.nVariables());

    std::iota(vec.begin(), vec.end(), 0);

    unsigned int i;

    auto str_tmp = BC_var.begin();

    PetscFunctionBegin;
    // going over all BC_ids
    for(i = 0; str_tmp != BC_var.end(); i++, str_tmp++)
    {
        std::vector<std::string> tmp = AssembleMassMatrix::split_string(*str_tmp, '-');

        // check if variable assigned in the input file exists for given simulation
        bool var_flg = 1;
        for(auto t = tmp.begin(); t != tmp.end(); ++t)
        {
            if(atoi(t->c_str()) >= _nl.nVariables())
                var_flg = 0;
        }

        // in case u havent put anything into input file, or u put too much
        if(*str_tmp == "-1" || var_flg == 0)
        {
            //std::cout<<"no_si"<<_nl.nVariables()<<std::endl;
            _dc_variables_id.push_back(vec);
        }
        else
        {
            unsigned int j;
            std::vector<int > one_BC_id;
            auto str_in = tmp.begin();
            for(j = 0; str_in != tmp.end(); j ++, str_in++)
            {
                one_BC_id.push_back(atoi(str_in->c_str()));
            }
            _dc_variables_id.push_back(one_BC_id);
        }
    }

    // check if u have same number of BC_ids in both parameters
    if(_dc_variables_id.size() != _dc_boundary_id.size())
    {
        _dc_variables_id.clear();
        for(auto i = 0; i != _dc_boundary_id.size(); i++)
        {
            _dc_variables_id.push_back(vec);
        }
    }

    // print out what is considered for zero-ing
    std::cout<<" ------ BC CONDITIONS Begin ------ \n";
    unsigned int t = 0;
    //std::cout<<"_dc_variables_id.begin()"<<_dc_variables_id.size()<<std::endl;
    for(auto i = _dc_variables_id.begin(); i != _dc_variables_id.end();  t++, i++)
    {
        std::cout<<"\n BC_id:  "<< _dc_boundary_id[t] << "   var_ids:  ";
        std::for_each(i->begin(), i->end(), [](int i){ std::cout << i << "  " ; });
    }
    std::cout<<" ------ BC CONDITIONS End ------ \n";
}


    

     
std::vector<std::string>
AssembleMassMatrix::split_string(const std::string & s, char delim)
{

      std::vector<std::string> v;

      if (s.length() == 0)
        std::cerr << "Got an empty string. Split_string(...) is confused. \n";

      auto i = 0;
      auto pos = s.find(delim);
      while (pos != std::string::npos)
      {
        v.push_back(s.substr(i, pos - i));
        i = ++pos;
        pos = s.find(delim, pos);

        if (pos == std::string::npos)
          v.push_back(s.substr(i, s.length()));
      }

      if (v.size() == 0) // if only one word is in the string
        v.push_back(s);

      return v;
}


