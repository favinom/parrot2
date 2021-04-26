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

  void getPermeability(std::vector<Point> const & p , std::vector<Real> & permeability);

  void getPermeabilityPoint(Point const & p , Real & permeability);

  void getPorosity(std::vector<Point> const & p , std::vector<Real> & porosity);

  void getPorosityPoint(Point const & p , Real & porosity);

protected:
    enum OperationType
    {
        FIRST,
        MAX,
        MIN,
        HARMONIC,
        ARITHMETIC
    };

    virtual void computeQpProperties();

    MooseEnum         const _operation_type;

    MeshGeneratorName       _meshGeneratorName;
    bool              const _hasMeshGenerator;

    Real              const _permeabilityBackInput;
    std::vector<Real>       _permeabilityFracInput;

    Real              const _porosityBackInput;
    Real              const _porosityFracInput;

    VariableGradient  const &_gradP;

    //Real const _poroInput;
    //Real _condFracture;
    
    //Real _dim;
    //RealTensorValue _id;
    //RealVectorValue _u_elem;

    //bool const _isPressureValid;
    //bool const _conservativeScheme;

    MaterialProperty<Real> &  _Kscalar;
    MaterialProperty<Real> &  _phi;
    MaterialProperty<RealVectorValue> &_U;

};
