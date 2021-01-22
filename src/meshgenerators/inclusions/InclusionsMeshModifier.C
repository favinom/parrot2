//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InclusionsMeshModifier.h"

#include <set>

#include "CastUniquePointer.h"

registerMooseObject("parrot2App", InclusionsMeshModifier);

defineLegacyParams(InclusionsMeshModifier);

InputParameters InclusionsMeshModifier::validParams()
{
	InputParameters params = MeshGenerator::validParams();

	params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");

	params.addRequiredParam<int>("fn", "number of fractures");
  
	params.addParam<int>("nx_min", 0, "");
	params.addParam<int>("nx_max", 0, "");
	params.addParam<int>("ny_min", 0, "");
	params.addParam<int>("ny_max", 0, "");
	params.addParam<int>("nz_min", 0, "");
	params.addParam<int>("nz_max", 0, "");
  
	params.addParam<Real>("Lx", 0.0, "Lx");
	params.addParam<Real>("Ly", 0.0, "Ly");
	params.addParam<Real>("Lz", 0.0, "Lz");
  
	params.addRequiredParam< std::vector<Real> >("fx_string", "x-coordinates of center of fractures");
	params.addRequiredParam< std::vector<Real> >("fy_string", "y-coordinates of center of fractures");
	params.addParam        < std::vector<Real> >("fz_string", "z-coordinates of center of fractures");
	params.addRequiredParam< std::vector<Real> >("fa1_string", "rotation along z axis");
	params.addParam        < std::vector<Real> >("fa2_string", "rotation along y axis");
	params.addParam        < std::vector<Real> >("fa3_string", "rotation along x axis");
	params.addRequiredParam< std::vector<Real> >("fd1_string", "fracture dimension 1");
	params.addRequiredParam< std::vector<Real> >("fd2_string", "fracture dimension 2");
	params.addParam        < std::vector<Real> >("fd3_string", "fracture dimension 3");

	return params;
}

InclusionsMeshModifier::InclusionsMeshModifier(const InputParameters & parameters) :
MeshGenerator(parameters),
_input(getMesh("input")),
_input_fn(getParam<int>("fn")),
_nx_min( getParam<int>("nx_min") ),
_nx_max( getParam<int>("nx_max") ),
_ny_min( getParam<int>("ny_min") ),
_ny_max( getParam<int>("ny_max") ),
_nz_min( getParam<int>("nz_min") ),
_nz_max( getParam<int>("nz_max") ),
_Lx( getParam<Real>("Lx") ),
_Ly( getParam<Real>("Ly") ),
_Lz( getParam<Real>("Lz") ),
_fx (getParam< std::vector<Real> >("fx_string")),
_fy (getParam< std::vector<Real> >("fy_string")),
_fa1(getParam< std::vector<Real> >("fa1_string")),
_fd1(getParam< std::vector<Real> >("fd1_string")),
_fd2(getParam< std::vector<Real> >("fd2_string"))
{

	if (parameters.isParamValid("fz_string")  &&
	 	parameters.isParamValid("fa2_string") &&
	 	parameters.isParamValid("fa3_string") &&
	 	parameters.isParamValid("fd3_string") )
	 {
        
	 	_fz =getParam< std::vector<Real> >("fz_string");
	 	_fa2=getParam< std::vector<Real> >("fa2_string");
	 	_fa3=getParam< std::vector<Real> >("fa3_string");
	 	_fd3=getParam< std::vector<Real> >("fd3_string");

	 	_console<<"Assuming dimension is 3\n";
	 	_fractureDim=3;
	 }
	 else if(!parameters.isParamValid("fz_string")  &&
	 		 !parameters.isParamValid("fa2_string") &&
	 		 !parameters.isParamValid("fa3_string") &&
	 		 !parameters.isParamValid("fd3_string") )
	 {
	 	_console<<"Assuming dimension is 2\n";
	 	_fractureDim=2;
	 }
	 else
	 {
	 	mooseError("Cannot detect the fracture dimension from input parameters!");
	 }

	 _total_fn=(_nz_max-_nz_min+1)*(_ny_max-_ny_min+1)*(_nx_max-_nx_min+1)*_input_fn;

	 if (_fx.size()  < _input_fn &&
	 	 _fy.size()  < _input_fn &&
	 	 _fa1.size() < _input_fn &&
	 	 _fd1.size() < _input_fn &&
	     _fd2.size() < _input_fn )
	 {
	 	mooseError("_fn is too large\n");
	 }
	 _fx.resize(_input_fn);
	 _fy.resize(_input_fn);
	 _fa1.resize(_input_fn);
	 _fd1.resize(_input_fn);
 	 _fd2.resize(_input_fn);

	 if ( _fractureDim==3 )
	 {
	 	if( _fz .size() < _input_fn &&
	 	    _fa2.size() < _input_fn &&
	 	    _fa3.size() < _input_fn &&
	 	    _fd3.size() < _input_fn )
	 	{
	 		mooseError("_fn is too large in 3D\n");
	 	}

	 	_fz .resize(_input_fn);
	 	_fa2.resize(_input_fn);
	 	_fa3.resize(_input_fn);
	 	_fd3.resize(_input_fn);
	 }

	std::vector< RealVectorValue > center;
	center.resize(_input_fn);
	std::vector< RealVectorValue > rotation;
	rotation.resize(_input_fn);
	std::vector< RealVectorValue > dimension;
	dimension.resize(_input_fn);

	for (int i=0; i<_input_fn; ++i)
	{
		center.at(i)(0)=_fx.at(i);
		center.at(i)(1)=_fy.at(i);
		center.at(i)(2)=0.0;

		rotation.at(i)(0)=_fa1.at(i)/180.0*pi;
		rotation.at(i)(1)=0.0;
		rotation.at(i)(2)=0.0;

		dimension.at(i)(0)=_fd1.at(i);
		dimension.at(i)(1)=_fd2.at(i);
		dimension.at(i)(2)=0.0;

		if (_fractureDim==3)
		{
			center.at(i)(2)   =_fz .at(i);
			rotation.at(1)    =_fa2.at(i)/180.0*pi;
			rotation.at(2)    =_fa3.at(i)/180.0*pi;
			dimension.at(i)(2)=_fd3.at(i);
		}
	
	}

	_boxInclusions.resize(_total_fn);
	int counter=0;
	
	for (int k=_nz_min; k<=_nz_max; ++k)
	{
	 	for (int j=_ny_min; j<=_ny_max; ++j)
	 	{
	 		for (int i=_nx_min; i<=_nx_max; ++i)
	 		{
	 			RealVectorValue shift(_Lx*(1.0*i),_Ly*(1.0*j),_Lz*(1.0*k));
	 			for (int f=0; f<_input_fn; ++f)
	 			{
	 				_boxInclusions.at(counter).Initialize(center.at(f)+shift,rotation.at(f),dimension.at(f),_fractureDim);
	 				++counter;
	 			}
	 		}
	 	}
	}
	
	/*
	std::cout<<_boxInclusions.size()<<std::endl;
	for (int i=0; i<_boxInclusions.size() ; ++i)
	{
	  	_boxInclusions.at(i).Print();
	}

	*/
}


std::unique_ptr<MeshBase> InclusionsMeshModifier::generate()
{
	std::unique_ptr<MeshBase> mesh = std::move(_input);
	return dynamic_pointer_cast<MeshBase>(mesh);
}

bool InclusionsMeshModifier::isInside(RealVectorValue const & point, Real bound) const
{
    for (int i=0; i<_total_fn; ++i)
    {
        if ( _boxInclusions.at(i).IsInside(point,bound) )
            return true;

    }
    return false;
}

std::vector<unsigned int> InclusionsMeshModifier::whichIsInside(RealVectorValue const & point, Real bound) const
{
	std::vector<unsigned int> output;

    for (int i=0; i<_total_fn; ++i)
    {
        if ( _boxInclusions.at(i).IsInside(point,bound) )
            output.push_back(i);

    }
    return output;
}


bool InclusionsMeshModifier::DoesIntersect(Elem const & elem) const
{
    for (int i=0; i<_total_fn; ++i)
    {
        if ( _boxInclusions.at(i).DoesIntersect(elem) )
            return true;
    }
    return false;
}

bool InclusionsMeshModifier::IsOnBoundary (Elem const & elem) const
{
	std::set<int> boundarySet;
	std::set<int> interiorSet;
	
    for (int i=0; i<_total_fn; ++i)
    {
        if (_boxInclusions.at(i).IsOnBoundary(elem) )
		{
			boundarySet.insert(i);
		}
    }
	
    for (int i=0; i<_total_fn; ++i)
    {
		if ( _boxInclusions.at(i).IsInside(elem) )
		{
			if ( boundarySet.count(i)==0 )
			{
				interiorSet.insert(i);
			}
		}
	}
	
	if ( !boundarySet.empty() && interiorSet.empty())
		return true;
	else
		return false;

}
