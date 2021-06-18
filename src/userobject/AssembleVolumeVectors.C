//* This file is part of the MOOSE framework
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AssembleVolumeVectors.h"
#include "InclusionsMeshModifier.h"




registerMooseObject("parrot2App", AssembleVolumeVectors);


template <>
InputParameters
validParams<AssembleVolumeVectors>()
{
  InputParameters params = validParams<GeneralUserObject>();
  
  params.addParam<MeshGeneratorName>("inclusions_list","inclusions_list");
  return params;
}

AssembleVolumeVectors::AssembleVolumeVectors(const InputParameters & parameters) :
GeneralUserObject(parameters),
   _hasMeshGenerator(parameters.isParamValid("inclusions_list"))
{

  if (_hasMeshGenerator)
  {
    _meshGeneratorName    =getParam<MeshGeneratorName>("inclusions_list");
  }

}

void AssembleVolumeVectors::execute()
{
  _console << "AssembleVolumeVectors::execute() begin"  << std::endl;

  //MeshModifier     const * _myMeshModifier_ptr;
  //RegionUserObject const * _myRegionUserObject_ptr;

  // if (_hasMeshModifier)
  // {
  //   MeshModifier const & _myMeshModifier( _app.getMeshModifier( _meshModifierName.c_str()) );
  //   //_myMeshModifier_ptr=&_myMeshModifier;
  //   RegionUserObject const & _regionUserObject( dynamic_cast<RegionUserObject const &>(_myMeshModifier) );
  //   _myRegionUserObject_ptr=&_regionUserObject;
  // }

  MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
  
  InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );
  


  DofMap   const & dof_map = _fe_problem.getNonlinearSystemBase().dofMap();
  auto &comm = _fe_problem.es().get_mesh().comm();

  _fn =   inclusionsMeshModifier.get_total_fn();
  volumes.resize(_fn);

  for (int i=0; i<_fn; ++i)
  {
    volumes[i]=new PetscVector<Number>(comm);
    volumes[i][0].init( dof_map.n_dofs() , dof_map.n_local_dofs() );
  }

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

  QBase const * const & qbase(_assembly.qRule());

  std::unique_ptr<QBase> qrule( QBase::build (qbase->type(),dim,qbase->get_order()));

  fe->attach_quadrature_rule (qrule.get());
  
  const std::vector<Real>& JxW      = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  //const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  const std::vector<Point>& q_points = fe->get_xyz();

  std::vector<dof_id_type> dof_indices;
  // this has to be resized
  std::vector< DenseVector<Number> > Re;
  std::vector< std::vector<Number> > coeff;

  Re.resize(_fn);
  coeff.resize(_fn);
  //coeff.zero();

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  MeshBase::const_element_iterator const end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem * elem = *el;
    fe->reinit (elem);
    dof_map.dof_indices(elem, dof_indices);

    int const loc_n=dof_indices.size();
    for (int i=0;i<_fn;++i)
    {
      Re.at(i).resize(loc_n);
      Re.at(i).zero();
      coeff.at(i).resize( qrule->n_points() );
    }

    Real ze=0.0;
    for (int qp=0; qp<qrule->n_points(); ++qp)
    {
        
      std::vector<unsigned int> which = inclusionsMeshModifier.whichIsInside( q_points[qp] );   
      // for (int i=0; i<_fn; ++i){
        for (int j=0; j<which.size(); ++j){     
          int d_id = which.at(j);
          // if(d_id == i){
            coeff.at(d_id).at(qp) = 1.0;
          // }
          // else{
          //   coeff.at(i).at(qp) = 0.0;
          // }
        //}
      }
    }

    for (unsigned int i=0; i<phi.size(); i++)
    {
      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
      {
        for (unsigned int j=0; j<_fn; ++j )
          Re.at(j)(i) += coeff.at(j).at(qp) * JxW[qp] * phi[i][qp];
      }
    }

    for (int i=0; i<_fn; ++i)
    {
      volumes[i][0].add_vector(Re.at(i), dof_indices);
    }
  }

  for (int i=0; i<_fn; ++i)
  {
    volumes[i][0].close();
    //std::cout<<volumes.at(i)[0].sum()<<std::endl<<std::endl;
  }

}

AssembleVolumeVectors::~AssembleVolumeVectors()
{
  
  MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
  
  InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );

  int fn =  inclusionsMeshModifier.get_fn();

  for (int i=0; i<fn; ++i)
    delete volumes[i];
}
