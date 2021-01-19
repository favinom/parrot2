//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InclusionRefinement.h"

#include "CastUniquePointer.h"

#include "libmesh/implicit_system.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/dof_map.h"

registerMooseObject("parrot2App", InclusionRefinement);

defineLegacyParams(InclusionRefinement);

InputParameters InclusionRefinement::validParams()
{
    InputParameters params = MeshGenerator::validParams();

    params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");

	params.addRequiredParam<MeshGeneratorName>("inclusions_list","inclusions_list");
	params.addParam< bool >("doBoundaryRefinement","doBoundaryRefinement");
	params.addParam< std::vector<int> >("refinements","refinements");
	params.addParam<std::string>("outputFileName","outputFileName");
    
    
	return params;
}

InclusionRefinement::InclusionRefinement(const InputParameters & parameters) :
MeshGenerator(parameters),
_input(getMesh("input")),
_inclusionsListName(getParam<MeshGeneratorName>("inclusions_list")),
_myMeshGenerator(_app.getMeshGenerator(_inclusionsListName)),
_inclusionsList( dynamic_cast<InclusionsMeshModifier const &>(_myMeshGenerator) ),
_hasRefinementVector( isParamValid("refinements") ),
_hasOutputFileName( isParamValid("outputFileName") )
{
	if (_hasRefinementVector)
	{
		_refinements=getParam<std::vector<int> >("refinements");
		if ( _refinements.size()%2 != 0 )
		{
			_console<<"The size of refinements has to be even, exiting...\n";
			exit(1);
		}
	}
	else
		_console<<"Maybe print mesh\n";

	if (_hasOutputFileName)
	{
		_outputFileName=getParam<std::string >("outputFileName");
	}
	if ( isParamValid("doBoundaryRefinement") )
	{
		_doBoundary=getParam<bool>("doBoundaryRefinement");
	}
	else
		_doBoundary=false;	
}

std::unique_ptr<MeshBase> InclusionRefinement::generate()
{
	std::cout<<"InclusionRefinement::generate()"<<std::endl;
	std::unique_ptr<MeshBase> mesh = std::move(_input);

	_meshBase=mesh.get();
	_unstructuredMesh=dynamic_cast<UnstructuredMesh *>(_meshBase);

	//TEST FOR PERIODIC // ONLY 2D	
	Real xmin(1e9),xmax(-1e9),ymin(1e9),ymax(-1e9);

	MeshBase::const_node_iterator     nd=_unstructuredMesh->local_nodes_begin();
	MeshBase::const_node_iterator end_nd=_unstructuredMesh->local_nodes_end  ();
	for ( ; nd != end_nd ; ++nd)
	{
		Node * node = *nd;
		xmin=std::min(xmin,node[0](0));
		ymin=std::min(ymin,node[0](1));
		xmax=std::max(xmax,node[0](0));
		ymax=std::max(ymax,node[0](1));
	}
	std::cout<<"min=("<<xmin<<","<<ymin<<")"<<std::endl;
	std::cout<<"max=("<<xmax<<","<<ymax<<")"<<std::endl;
	
	EquationSystems equation_systems (_unstructuredMesh[0]);
	ImplicitSystem & system = equation_systems.add_system<ImplicitSystem> ("system");
	system.add_variable ("variable", FIRST);
	
	DofMap & dof_map = system.get_dof_map();
	
	PeriodicBoundary horz(RealVectorValue(xmax-xmin, 0., 0.));
	horz.myboundary = 3;
	horz.pairedboundary = 1;
	dof_map.add_periodic_boundary(horz);
	
	PeriodicBoundary vert(RealVectorValue(0., ymax-ymin, 0.));
	vert.myboundary = 0;
	vert.pairedboundary = 2;
	dof_map.add_periodic_boundary(vert);
	
	equation_systems.init ();	
	//END TEST
	meshRefinement=new MeshRefinement(_unstructuredMesh[0]);
	meshRefinement[0].set_periodic_boundaries_ptr(dof_map.get_periodic_boundaries());

	if (_hasRefinementVector)
		doRefine(_refinements);

	if (_hasOutputFileName)
	{
		_console<<"\nWriting mesh file\n";
		_unstructuredMesh->write(_outputFileName);
		_console<<"Done!\n";
	}

	std::cout<<"InclusionRefinement::generate()"<<std::endl;

	return dynamic_pointer_cast<MeshBase>(mesh);
}

void InclusionRefinement::doAMR()
{
	MeshBase::const_element_iterator     el=_unstructuredMesh->active_elements_begin();
	MeshBase::const_element_iterator end_el=_unstructuredMesh->active_elements_end  ();

	std::cout<<"  Looping over elements\n";
	for ( ; el != end_el ; ++el)
	{
		Elem * elem = *el;
		if (_doBoundary)
		{
			if ( _inclusionsList.IsOnBoundary(elem[0]) )
			{
				elem[0].set_refinement_flag(Elem::REFINE);
			}
		}
		else
		{
			if ( _inclusionsList.DoesIntersect(elem[0]) )
				elem[0].set_refinement_flag(Elem::REFINE);
		}
	}

	std::cout<<"  Actual refinement\n";
	meshRefinement->refine_elements();
	std::cout<<"  Done!\n";
}

void InclusionRefinement::doUMR(int i)
{
	meshRefinement->uniformly_refine(i);
}


void InclusionRefinement::doRefine(std::vector<int> const & refinements)
{
	std::cout<<"We are going to perform this sequence of refinements\n";
	for (int i=0; i<refinements.size(); ++++i)
	{
		std::cout<<refinements.at(i)<<" AMR and ";
		std::cout<<refinements.at(i+1)<<" UMR\n";
	}

	for (int i=0;i<refinements.size(); ++++i)
	{
		std::cout<<"\nAdaptivty loop " <<i/2<<": "<<refinements.at(i)<<" AMR and "<<refinements.at(i+1)<<" UMR"<<std::endl;
		for (int j=0; j<refinements.at(i); ++j)
		{
			std::cout<<" Doing AMR "<<j+1<<std::endl;
			doAMR();
		}

		std::cout<<" Doing "<<refinements.at(i+1)<<" UMR"<<std::endl;
		if (refinements.at(i+1)>0)
			doUMR(refinements.at(i+1));
	}

}
