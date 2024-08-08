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

  void getPermeabilityId(std::vector<Point> const & p , std::vector<Real> & permeability, int & ID);

  void getPermeabilityPoint(Point const & p , Real & permeability);

  void getPermeabilityPointId(Point const & p , Real & permeability, int & ID);

  void getPorosity(std::vector<Point> const & p , std::vector<Real> & porosity);

  void getPorosityId(std::vector<Point> const & p , std::vector<Real> & porosity, int & ID);

  void getPorosityPoint(Point const & p , Real & porosity);

  void getPorosityPointId(Point const & p , Real & porosity, int & ID);

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

    std::vector<int> _vector_p;
    std::vector<Real> _vector_value_k;
    std::vector<Real> _vector_value_phi;
    bool _multisubdomain;

};
