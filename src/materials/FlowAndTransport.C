#include "FlowAndTransport.h"
#include "InclusionsMeshModifier.h"

registerMooseObject("parrot2App", FlowAndTransport);

defineLegacyParams(FlowAndTransport);

InputParameters FlowAndTransport::validParams()
{
	InputParameters params = Material::validParams();

	params.addParam<MeshGeneratorName>("inclusions_list","inclusions_list");
	
	params.addRequiredParam<Real>("k","k");
	params.addParam<Real>("kFrac","kFrac");

	//params.addRequiredParam<bool>("conservative","use a conservative scheme?");
	//params.addCoupledVar("pressure","The gradient of this variable will be used as the velocity vector.");

	//params.addRequiredParam<Real>("phi","phi");
	//params.addParam<Real>("phiFrac","phiFrac");

	return params;
}

FlowAndTransport::FlowAndTransport(const InputParameters &parameters) :
Material(parameters),
_hasMeshGenerator(parameters.isParamValid("inclusions_list")),
_permeabilityBackInput(getParam<Real>("k")),
_Kscalar(declareProperty<Real>("conductivityProperty"))/*,
_poroInput(getParam<Real>("phi")),
_dim(_mesh.dimension()),
_isPressureValid(parameters.isParamValid("pressure")),
_conservativeScheme(getParam<bool>("conservative")),
_gradP(_isPressureValid ? coupledGradient("pressure"): _grad_zero),
_poro(declareProperty<Real>("Porosity")),
_K(declareProperty<RealTensorValue>("PermeabilityTensor")),

_U(declareProperty<RealVectorValue>("VelocityVector"))*/
{
	// _id=RealTensorValue(0.0,0.0,0.0,
	// 	0.0,0.0,0.0,
	// 	0.0,0.0,0.0);

	// for (int i=0; i<_dim; ++i)
	// 	_id(i,i)=1.0;

	if ( _hasMeshGenerator && !isParamValid("kFrac"))
	{
			mooseError("At least one of the following is missing: inclusions_list, kFrac");
	}

	if (_hasMeshGenerator)
	{
		_meshGeneratorName    =getParam<MeshGeneratorName>("inclusions_list");
  		_permeabilityFracInput=getParam<Real>("kFrac");
  		//_poroFracture=getParam<Real>("phiFrac");
	}

}

void
FlowAndTransport::computeQpProperties()
{
	_Kscalar[_qp]=_permeabilityBackInput;
	if (_hasMeshGenerator)
	{
	 	MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
	 	InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );
	 	if (inclusionsMeshModifier.isInside(_q_point[_qp]))
	 	{
	// 		_poro[_qp]   =_poroFracture;
	 		_Kscalar[_qp]=_permeabilityFracInput;
	// 		_K[_qp]      =_condFracture * _id;
	 	}
	}

	// _poro[_qp]   =_poroInput;
	
	// _K[_qp]      =_condInput * _id;


	// if (_isPressureValid)
	// {
	// 	_U[_qp] =  -1.0 * _K[_qp] * _gradP[_qp];
	// }
	// else
	// {
	// 	for (int i=0; i<3; ++i)
	// 		_U[_qp](i)=0.0/0.0;
	// }

	// if (_conservativeScheme && _qp==_qrule->n_points()-1)
	// {
	// 	RealVectorValue b;
	// 	RealTensorValue A;
	// 	b.zero();
	// 	A.zero();
	// 	for (int i=0; i<_dim; ++i)
	// 	{
	// 		for (int qp=0; qp<_qrule->n_points(); ++qp)
	// 		{
	// 			b(i)+=_JxW[qp]*_gradP[qp](i);
	// 			A(i,i)+=-1.0*_JxW[qp]/_Kscalar[qp];
	// 		}
	// 		_u_elem(i)=b(i)/A(i,i);
	// 	}
	// 	for (int qp=0; qp<_qrule->n_points(); ++qp)
	// 	{
	// 		_U[qp]=_u_elem;
	// 	}
	// }
}
