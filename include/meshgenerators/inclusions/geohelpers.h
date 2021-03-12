#include "libmesh/vector_value.h"

using namespace libMesh;

bool doesEdgeIntersectElement_2D(RealVectorValue const & p1, RealVectorValue const & p2, RealVectorValue const & hmin, RealVectorValue const & hmax);
Real CalcY_2D(Real const & xval, RealVectorValue const & p1, RealVectorValue const & p2);
Real CalcX_2D(Real const & yval, RealVectorValue const & p1, RealVectorValue const & p2);


bool doesEdgeIntersectElement_3D(RealVectorValue const & p1, RealVectorValue const & p2, RealVectorValue const & hmin, RealVectorValue const & hmax); 

Real CalcXY_3D(Real const & zval, RealVectorValue const & p1, RealVectorValue const & p2, int const dir);
Real CalcXZ_3D(Real const & yval, RealVectorValue const & p1, RealVectorValue const & p2, int const dir);
Real CalcYZ_3D(Real const & xval, RealVectorValue const & p1, RealVectorValue const & p2, int const dir);

//void CalcXY_3D(Real const & zval, RealVectorValue const & p1, RealVectorValue const & p2, Real & xval, Real & yval);
