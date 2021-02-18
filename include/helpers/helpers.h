#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"

using namespace libMesh;

void copyVariable(EquationSystems const & eqIn, std::string const variableIn, EquationSystems & eqOut, std::string const variableOut);