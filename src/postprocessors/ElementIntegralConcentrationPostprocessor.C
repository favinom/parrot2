//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegralConcentrationPostprocessor.h"
#include "InclusionsMeshModifier.h"

registerMooseObject("parrot2App", ElementIntegralConcentrationPostprocessor);

template <>
InputParameters
validParams<ElementIntegralConcentrationPostprocessor>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addParam<MeshGeneratorName>("inclusions_list","inclusions_list");
  params.addRequiredParam<unsigned int>("fractureRegionId","fractureRegionId");
  return params;
}

ElementIntegralConcentrationPostprocessor::ElementIntegralConcentrationPostprocessor(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
   _hasMeshGenerator(parameters.isParamValid("inclusions_list")), 
   _regionId(getParam<unsigned int>("fractureRegionId")),
   _u(coupledValue("variable"))
{
  addMooseVariableDependency(mooseVariable());
  if (_hasMeshGenerator)
  {
    _meshGeneratorName    =getParam<MeshGeneratorName>("inclusions_list");
  }

}


Real
ElementIntegralConcentrationPostprocessor::computeQpIntegral()
{
  // params.addRequiredParam<std::string>("fractureMeshModifier","fractureMeshModifier");
  // _meshModifierName(getParam<std::string>("fractureMeshModifier")),
    MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
    InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );

    std::vector<unsigned int> which = inclusionsMeshModifier.whichIsInside( _q_point[_qp] );

    unsigned int check = 0;

    bool bool_rg = false;

    int i = 0;

    while (i<which.size() && bool_rg==false)
    {
      
      check = which.at(i);
      //std::cout<<"check="<<check<<std::endl;

      if(check ==_regionId){
        bool_rg=true;
      }
      
      i++;


     }

    //std::cout<<bool_rg<<std::endl;   

    if (bool_rg==true)
    return _u[_qp];
    else
    return 0.0;
}

