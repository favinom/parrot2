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
#include "MeshGenerator.h"
#include "FEProblem.h"
//#include "MooseEnum.h"

#include "libmesh/distributed_mesh.h"

#include "InclusionsMeshModifier.h"


// Forward declerations
class InclusionRefinement;

template <>
InputParameters validParams<InclusionRefinement>();

class InclusionRefinement : public MeshGenerator
{
public:

	static InputParameters validParams();

    InclusionRefinement(const InputParameters & parameters);
   	
   	virtual std::unique_ptr<MeshBase> generate() override;
    void doAMR();
    void doUMR(int i);
    void doRefine(std::vector<int> const & refinements);
    
protected:

	std::unique_ptr<MeshBase> & _input;

	MeshGeneratorName      const   _inclusionsListName;
	MeshGenerator          const & _myMeshGenerator;
	InclusionsMeshModifier const & _inclusionsList;

	bool _hasRefinementVector;
	std::vector<int> _refinements;

	bool _hasOutputFileName;
	std::string _outputFileName;

	bool _doBoundary;

	MeshBase * _meshBase;
	UnstructuredMesh * _unstructuredMesh;
	MeshRefinement * meshRefinement;
	
};
