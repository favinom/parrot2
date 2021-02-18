#include "helpers.h"

void copyVariable(EquationSystems const & eqIn, std::string const variableIn, EquationSystems & eqOut, std::string const variableOut)
{
	int id_in;
	int id_out;

	DofMap const * df_in;
	DofMap const * df_out;

	NumericVector<Number> * sol_in;
	NumericVector<Number> * sol_out;

	bool foundIn=false;
	bool foundOut=false;

	MeshBase const & meshBase( eqOut.get_mesh() );

	for (int i=0; i<eqIn.n_systems(); ++i)
	{
		System const & es=eqIn.get_system<System>(i);
		if ( es.has_variable(variableIn) )
		{
			id_in=es.variable_number(variableIn);
			DofMap const & df_in_ref=es.get_dof_map();
			df_in=&df_in_ref;
			std::unique_ptr<NumericVector<Number>> const & solutionPointer(es.solution);
			sol_in=solutionPointer.get();
			std::cout<<std::endl<<std::endl<<"In the equation systems of input,\n";
			std::cout<<"found a system with name "<<es.name()<<std::endl;
			std::cout<<"which contains a variable with the id "<<id_in<<std::endl<<std::endl;
			foundIn=true;
			break;
		}
	}

	for (int i=0; i<eqOut.n_systems(); ++i)
	{
		System const & es=eqOut.get_system<System>(i);
		if ( es.has_variable(variableOut) )
		{
			id_out=es.variable_number(variableOut);
			DofMap const & df_out_ref=es.get_dof_map();
			df_out=&df_out_ref;
			std::unique_ptr<NumericVector<Number>> const & solutionPointer(es.solution);
			sol_out=solutionPointer.get();
			std::cout<<std::endl<<std::endl<<"In the equation systems of output,\n";
			std::cout<<"found a system with name "<<es.name()<<std::endl;
			std::cout<<"which contains a variable with the id "<<id_in<<std::endl<<std::endl;
			foundOut=true;
			break;
		}
	}

	if (foundOut==false || foundIn==false)
	{
		std::cout<<"I am not proceeding to copy, as input or output variable has not been found.";
	}
	else
	{
		MeshBase::const_node_iterator       node_it     = meshBase.local_nodes_begin();
		MeshBase::const_node_iterator const node_it_end = meshBase.local_nodes_end();

		std::cout<<sol_in[0].size()<<" "<<sol_in[0].initialized()<<std::endl;
		std::cout<<sol_out[0].size()<<" "<<sol_out[0].initialized()<<std::endl;

		for ( ; node_it != node_it_end; ++node_it)
    	{
        	const Node * node = *node_it;
        	std::vector<dof_id_type> dof_in;
        	std::vector<dof_id_type> dof_out;
        	
        	df_in [0].dof_indices (node, dof_in ,id_in );
        	df_out[0].dof_indices (node, dof_out,id_out);

        	if (dof_in.size()!=1 && dof_out.size()!=1 )
        	{
        		std::cout<<"Error: dof_in.size()!=1 && dof_out.size()!=1, exiting...\n";
        		exit(1);
        	}
        	Number value;
        	sol_in [0].get(dof_in       ,&value);
        	sol_out[0].set(dof_out.at(0), value);
        }
        sol_out[0].close();
        eqOut.update();
    }
}
