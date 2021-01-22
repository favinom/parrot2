#pragma once

#include "Material.h"

class FlowAndTransport;

template <>
InputParameters validParams<FlowAndTransport>();

class FlowAndTransport : public Material
{
public:
  static InputParameters validParams();

  FlowAndTransport(const InputParameters &parameters);

protected:
  virtual void computeQpProperties();


protected:

    MeshGeneratorName       _meshGeneratorName;
    bool              const _hasMeshGenerator;

    Real const _permeabilityBackInput;
    Real       _permeabilityFracInput;

    //Real const _poroInput;
    //Real _condFracture;
    
    //Real _dim;
    //RealTensorValue _id;
    //RealVectorValue _u_elem;

    //bool const _isPressureValid;
    //bool const _conservativeScheme;

    //const VariableGradient &_gradP;
  
    //MaterialProperty<Real> &_poro;
    //MaterialProperty<RealTensorValue> &_K;
    MaterialProperty<Real> &_Kscalar;
    //MaterialProperty<RealVectorValue> &_U;

};
