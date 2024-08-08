#include "FlowAndTransport.h"
#include "InclusionsMeshModifier.h"

registerMooseObject("parrot2App", FlowAndTransport);

defineLegacyParams(FlowAndTransport);

InputParameters FlowAndTransport::validParams()
{
	InputParameters params = Material::validParams();

	params.addParam<MeshGeneratorName>("inclusions_list","inclusions_list");
	
	// This Enum contains all possible operations possible when different fractures intersect
	MooseEnum operation_type("first max min harmonic arithmetic","first");
	params.addParam<MooseEnum>("operation_type", operation_type, "Which operation?");

	params.addRequiredParam<Real>("k","k");
	params.addParam< std::vector<Real> >("kFrac","kFrac");

	params.addCoupledVar("pressure","The gradient of this variable will be used as the velocity vector.");

	params.addParam<Real>("phi",0.0,"phi");
	params.addParam<Real>("phiFrac",0.0,"phiFrac");
	params.addParam<std::vector<int>>("block_id","block_id");
	params.addParam<std::vector<Real>>("value_phi","value_phi");
	params.addParam<std::vector<Real>>("value_K","value_K");
	params.addParam<bool>("multisubdomain",false,"multisubdomain");


	return params;
}

FlowAndTransport::FlowAndTransport(const InputParameters &parameters) :
Material(parameters),
_operation_type(getParam<MooseEnum>("operation_type")),
_hasMeshGenerator(parameters.isParamValid("inclusions_list")),
_permeabilityBackInput(getParam<Real>("k")),
_porosityBackInput(getParam<Real>("phi")),
_porosityFracInput(getParam<Real>("phiFrac")),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"): _grad_zero),
_Kscalar(declareProperty<Real>("conductivityProperty")),
_phi(declareProperty<Real>("porosityProperty")),
_U(declareProperty<RealVectorValue>("velocityProperty")),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value_k(getParam<std::vector<Real>>("value_K")),
_vector_value_phi(getParam<std::vector<Real>>("value_phi")),
_multisubdomain(getParam<bool>("multisubdomain"))
{

	if ( _hasMeshGenerator && !isParamValid("kFrac"))
	{
			mooseError("At least one of the following is missing: inclusions_list, kFrac");
	}

	if (_hasMeshGenerator)
	{
		_meshGeneratorName    =getParam<MeshGeneratorName>("inclusions_list");
  		_permeabilityFracInput=getParam< std::vector<Real> >("kFrac");
  		
  		MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
	 	InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );

	 	if ( _permeabilityFracInput.size()>inclusionsMeshModifier.get_fn() )
	 		_permeabilityFracInput.resize( inclusionsMeshModifier.get_fn() );

	 	if ( _permeabilityFracInput.size()==1 )
	 	{
	 		_permeabilityFracInput.resize( inclusionsMeshModifier.get_fn() , _permeabilityFracInput.at(0) );
	 	}

	 	if ( _permeabilityFracInput.size()!=inclusionsMeshModifier.get_fn() )
	 	{
	 		mooseError("_permeabilityFracInput.size()!=inclusionsMeshModifier.get_fn()");
	 	}

	 	int fn=inclusionsMeshModifier.get_fn();
	 	int tfn=inclusionsMeshModifier.get_total_fn();
	 	for (int i=fn; i<tfn; ++i)
	 	{
	 		int  index=i%fn;
	 		Real value=_permeabilityFracInput.at(index);
	 		_permeabilityFracInput.push_back(value);
	 	}
	}

}

void
FlowAndTransport::computeQpProperties()
{
	getPermeabilityPoint( _q_point[_qp] , _Kscalar[_qp] );
	getPorosityPoint    ( _q_point[_qp] , _phi[_qp]     );
	_U[_qp]=-1.0 * _Kscalar[_qp]*_gradP[_qp];
	//std::cout<<"FlowAndTransport::computeQpProperties()"<<	_gradP[_qp]<<std::endl;
}

void FlowAndTransport::getPermeability(std::vector<Point> const & p , std::vector<Real> & permeability)
{
	permeability.resize( p.size() );
	for (int i=0; i<p.size(); ++i)
	{
		getPermeabilityPoint( p.at(i) , permeability.at(i) );
	}
}



void FlowAndTransport::getPermeabilityId(std::vector<Point> const & p , std::vector<Real> & permeability, int & subdomain_id){

	permeability.resize( p.size() );
	
	for (int i=0; i<p.size(); ++i)
	{
		getPermeabilityPointId( p.at(i) , permeability.at(i), subdomain_id);
	}
}

void FlowAndTransport::getPorosity(std::vector<Point> const & p , std::vector<Real> & porosity)
{
	porosity.resize( p.size() );
	for (int i=0; i<p.size(); ++i)
	{
		getPorosityPoint( p.at(i) , porosity.at(i) );
	}
}

void FlowAndTransport::getPorosityId(std::vector<Point> const & p , std::vector<Real> & porosity, int & subdomain_id)
{
	porosity.resize( p.size() );
	for (int i=0; i<p.size(); ++i)
	{
		getPorosityPointId( p.at(i) , porosity.at(i), subdomain_id );
	}
}


void FlowAndTransport::getPermeabilityPoint(Point const & p , Real & permeability)
{


	permeability = 0.0;

	if (_multisubdomain){

		for(int ll=0; ll<_vector_p.size(); ll++)
		{
			if (_assembly.elem()->subdomain_id()==_vector_p[ll])
			{
				permeability = _vector_value_k[ll];

				//std::cout<<"domain_id: "<< _assembly.elem()->subdomain_id()<<"and "<<_assembly.currentSubdomainID() <<std::endl;
				
				//std::cout<<"permeability: "<<permeability<<std::endl;
			}
		}
	}
	else{
		permeability=_permeabilityBackInput;
	}

	if (_hasMeshGenerator)
	{
	  	MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
	  	InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );

	  	std::vector<unsigned int> which=inclusionsMeshModifier.whichIsInside( p );
	  	if ( which.size()>0 )
	 	{
	 		std::vector<Real> permLocal( which.size() );
	 		for (int i=0; i<which.size(); ++i)
	 		{
	 			permLocal.at(i)=_permeabilityFracInput.at( which.at(i) );
	 		}
	 		switch ( _operation_type )
	 		{
    			case FIRST:
    			{
    				permeability=permLocal.at(0);
    			}
        		break;
    			case MAX:
    			{
    				permeability=permLocal.at(0);
    				for (int i=1; i<permLocal.size(); ++i)
    				{
    					permeability=std::max(permeability,permLocal.at(i));
    				}
    			}
        		break;
        		case MIN:
    			{
    				permeability=permLocal.at(0);
    				for (int i=1; i<permLocal.size(); ++i)
    				{
    					permeability=std::min(permeability,permLocal.at(i));
    				}
    			}
        		break;
        		case HARMONIC:
    			{
    				Real result=0.0;
    				for (int i=0; i<permLocal.size(); ++i)
    					result+=1.0/permLocal.at(i);
    				permeability=permLocal.size()/result;
    				//if (permLocal.size()>1)
    				//	std::cout<<permeability<<std::endl<<std::endl;
    			}
        		break;
        		case ARITHMETIC:
    			{
    				Real result=std::accumulate(std::begin(permLocal), std::end(permLocal),0.0);
    				permeability=result/permLocal.size();
    			}
        		break;
    			default:
    			{
    				mooseError("no case");// code to be executed if n doesn't match any cases
    			}
			}// switch
	  	}
	}
}


void FlowAndTransport::getPermeabilityPointId(Point const & p , Real & permeability, int & subdomain_id)
{


	permeability = 0.0;

	//std::cout<<"domain_id_out: "<< subdomain_id <<std::endl;

	for(int ll=0; ll<_vector_p.size(); ll++)
	{
		if (subdomain_id ==_vector_p[ll])
		{
				permeability = _vector_value_k[ll];
				
				//std::cout<<"permeability: "<<permeability<<std::endl;
		}
	}
	
	

	if (_hasMeshGenerator)
	{
	  	MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
	  	InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );

	  	std::vector<unsigned int> which=inclusionsMeshModifier.whichIsInside( p );
	  	if ( which.size()>0 )
	 	{
	 		std::vector<Real> permLocal( which.size() );
	 		for (int i=0; i<which.size(); ++i)
	 		{
	 			permLocal.at(i)=_permeabilityFracInput.at( which.at(i) );
	 		}
	 		switch ( _operation_type )
	 		{
    			case FIRST:
    			{
    				permeability=permLocal.at(0);
    			}
        		break;
    			case MAX:
    			{
    				permeability=permLocal.at(0);
    				for (int i=1; i<permLocal.size(); ++i)
    				{
    					permeability=std::max(permeability,permLocal.at(i));
    				}
    			}
        		break;
        		case MIN:
    			{
    				permeability=permLocal.at(0);
    				for (int i=1; i<permLocal.size(); ++i)
    				{
    					permeability=std::min(permeability,permLocal.at(i));
    				}
    			}
        		break;
        		case HARMONIC:
    			{
    				Real result=0.0;
    				for (int i=0; i<permLocal.size(); ++i)
    					result+=1.0/permLocal.at(i);
    				permeability=permLocal.size()/result;
    				//if (permLocal.size()>1)
    				//	std::cout<<permeability<<std::endl<<std::endl;
    			}
        		break;
        		case ARITHMETIC:
    			{
    				Real result=std::accumulate(std::begin(permLocal), std::end(permLocal),0.0);
    				permeability=result/permLocal.size();
    			}
        		break;
    			default:
    			{
    				mooseError("no case");// code to be executed if n doesn't match any cases
    			}
			}// switch
	  	}
	}
}

void FlowAndTransport::getPorosityPoint(Point const & p , Real & porosity)
{
	porosity = 0.0;
	
	if (_multisubdomain){

		for(int ll=0; ll<_vector_p.size(); ll++){

			if (_assembly.elem()->subdomain_id()==_vector_p[ll])
			{
				porosity = _vector_value_phi[ll];

				//if(_assembly.elem()->subdomain_id()==1) std::cout<<"_assembly.elem()->subdomain_id()"<<_assembly.elem()->subdomain_id()<<" and porosity: "<<porosity<<std::endl;
			}
		}
	}
	else{
		porosity=_porosityBackInput;
	}

			
	if (_hasMeshGenerator)
	{
		MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
	  	InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );
	  	if (inclusionsMeshModifier.isInside( p ) )
	  	{
	  		porosity=_porosityFracInput;
	  	}
	}
}



void FlowAndTransport::getPorosityPointId(Point const & p , Real & porosity, int & subdomain_id)
{
	porosity = 0.0;

	for(int ll=0; ll<_vector_p.size(); ll++){
		if (subdomain_id ==_vector_p[ll])
		{
			porosity = _vector_value_phi[ll];
	    }
	}
	
			
	if (_hasMeshGenerator)
	{
		MeshGenerator          const & myMeshGenerator       ( _app.getMeshGenerator( _meshGeneratorName ) );
	  	InclusionsMeshModifier const & inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(myMeshGenerator) );
	  	if (inclusionsMeshModifier.isInside( p ) )
	  	{
	  		porosity=_porosityFracInput;
	  	}
	}
}