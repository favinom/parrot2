#include "geohelpers.h"

#define myeps 1e-14

bool doesEdgeIntersectElement_2D(RealVectorValue const & p1, RealVectorValue const & p2, RealVectorValue const & hmin, RealVectorValue const & hmax)
{

	Real const & left  =hmin(0);
	Real const & right =hmax(0);
	Real const & bottom=hmin(1);
	Real const & top   =hmax(1);
	
	Real max_p_x=std::max(p1(0),p2(0));
	Real min_p_x=std::min(p1(0),p2(0));

	Real max_p_y=std::max(p1(1),p2(1));
	Real min_p_y=std::min(p1(1),p2(1));

	if (right+myeps<min_p_x)
		return false;
	if (max_p_x+myeps<left)
		return false;
	if (top+myeps<min_p_y)
		return false;
	if (max_p_y+myeps<bottom)
		return false;

	Real xb=CalcX_2D(bottom,p1,p2);
	Real xt=CalcX_2D(top   ,p1,p2);
	Real yl=CalcY_2D(left  ,p1,p2);
	Real yr=CalcY_2D(right ,p1,p2);

	//int intersections=0;

	if( left-myeps<xb && xb<right+myeps )
	{
		//intersections=intersections+1;
		return true;
	}
	if( left-myeps<xt && xt<right+myeps )
	{
		//intersections=intersections+1;
		return true;
	}
	if( bottom-myeps<yl && yl<top+myeps )
	{
		//intersections=intersections+1;
		return true;
	}
	if( bottom-myeps<yr && yr<top+myeps )
	{
		//intersections=intersections+1;
		return true;
	}

	if ( hmin(0)-myeps<p1(0) && p1(0)<hmax(0)+myeps && hmin(1)-myeps<p1(1) && p1(1)<hmax(1)+myeps )
		return true;
	if ( hmin(0)-myeps<p2(0) && p2(0)<hmax(0)+myeps && hmin(1)-myeps<p2(1) && p2(1)<hmax(1)+myeps )
		return true;

	return false;
}

Real CalcY_2D(Real const & xval, RealVectorValue const & p1, RealVectorValue const & p2)
{
    Real x0=p1(0);
    Real x1=p2(0);
    Real y0=p1(1);
    Real y1=p2(1);
    if ( std::fabs(x0-x1)<myeps )
    	return NAN;

    return y0 + (xval - x0)*(y1 - y0)/(x1 - x0);
}

Real CalcX_2D(Real const & yval, RealVectorValue const & p1, RealVectorValue const & p2)
{
	Real x0=p1(0);
    Real x1=p2(0);
    Real y0=p1(1);
    Real y1=p2(1);

    if ( std::fabs(y0-y1)<myeps )
    	return NAN;

    return x0 + (yval - y0)*(x1 - x0)/(y1 - y0);
}


bool doesEdgeIntersectElement_3D(RealVectorValue const & p1, RealVectorValue const & p2, RealVectorValue const & hmin, RealVectorValue const & hmax)
{
	
	Real const left  =hmin(0);
	Real const right =hmax(0);
	Real const bottom=hmin(1);
	Real const top   =hmax(1);
	Real const front =hmin(2);
	Real const back  =hmax(2);

	Real const max_p_x=std::max(p1(0),p2(0));
	Real const min_p_x=std::min(p1(0),p2(0));

	Real const max_p_y=std::max(p1(1),p2(1));
	Real const min_p_y=std::min(p1(1),p2(1));

	Real const max_p_z=std::max(p1(2),p2(2));
	Real const min_p_z=std::min(p1(2),p2(2));

	if (right+myeps<min_p_x)
		return false;
	if (max_p_x+myeps<left)
		return false;
	if (top+myeps<min_p_y)
		return false;
	if (max_p_y+myeps<bottom)
		return false;
	if (back+myeps<min_p_z)
		return false;
	if (max_p_z+myeps<front)
		return false;
	
	Real xb=CalcXZ_3D(bottom,p1,p2,0);
	Real zb=CalcXZ_3D(bottom,p1,p2,2);
	Real xt=CalcXZ_3D(top   ,p1,p2,0);
	Real zt=CalcXZ_3D(top   ,p1,p2,2);
	
	
	Real yl=CalcYZ_3D(left  ,p1,p2,1);
	Real zl=CalcYZ_3D(left  ,p1,p2,2);
	Real yr=CalcYZ_3D(right ,p1,p2,1);
	Real zr=CalcYZ_3D(right ,p1,p2,2);

	Real xf=CalcXY_3D(front ,p1,p2,0);
	Real yf=CalcXY_3D(front ,p1,p2,1);
	Real x2=CalcXY_3D(back  ,p1,p2,0);
	Real y2=CalcXY_3D(back  ,p1,p2,1);

	//std::cout<<xb<<std::endl;
	if( left-myeps<xb && xb<right+myeps && front-myeps<zb && zb<back+myeps ) // 
		return true;
	if( left-myeps<xt && xt<right+myeps && front-myeps<zt && zt<back+myeps)
		return true;
	if( bottom-myeps<yl && yl<top+myeps && front-myeps<zl && zl<back+myeps)
		return true;
	if( bottom-myeps<yr && yr<top+myeps && front-myeps<zr && zr<back+myeps)
		return true;
	if( left-myeps<xf && xf<right+myeps && bottom-myeps<yf && yf<top+myeps)
		return true;
	if( left-myeps<x2 && x2<right+myeps && bottom-myeps<y2 && y2<top+myeps)
		return true;
	
	if ( hmin(0)-myeps<p1(0) && p1(0)<hmax(0)+myeps && hmin(1)-myeps<p1(1) && p1(1)<hmax(1)+myeps && hmin(2)-myeps<p1(2) && p1(2)<hmax(2)+myeps ) 
		return true;
	if ( hmin(0)-myeps<p2(0) && p2(0)<hmax(0)+myeps && hmin(1)-myeps<p2(1) && p2(1)<hmax(1)+myeps && hmin(2)-myeps<p2(2) && p2(2)<hmax(2)+myeps )
		return true;
	
	return false;

}


Real CalcXY_3D(Real const & zval, RealVectorValue const & p1, RealVectorValue const & p2, int const dir)
{
	Real const x0=p1(0);
	Real const y0=p1(1);
	Real const z0=p1(2);
	Real const x1=p2(0);
	Real const y1=p2(1);
	Real const z1=p2(2);

	if ( std::fabs(z0-z1)<myeps )
		return NAN;

	if (dir==0)
		return x0 + (zval - z0)/(z1 - z0)*(x1 - x0);
	if (dir==1)
		return y0 + (zval - z0)/(z1 - z0)*(y1 - y0);

	exit(1);
	return 0.0;
}

Real CalcXZ_3D(Real const & yval, RealVectorValue const & p1, RealVectorValue const & p2, int const dir)
{
	Real const x0=p1(0);
	Real const y0=p1(1);
	Real const z0=p1(2);
	Real const x1=p2(0);
	Real const y1=p2(1);
	Real const z1=p2(2);

	if ( std::fabs(y0-y1)<myeps )
		return NAN;

	if (dir==0)
		return x0 + (yval - y0)/(y1 - y0)*(x1 - x0);
	if (dir==2)
		return z0 + (yval - y0)/(y1 - y0)*(z1 - z0);

	exit(1);
	return 0.0;
}

Real CalcYZ_3D(Real const & xval, RealVectorValue const & p1, RealVectorValue const & p2, int const dir)
{
	Real const x0=p1(0);
	Real const y0=p1(1);
	Real const z0=p1(2);
	Real const x1=p2(0);
	Real const y1=p2(1);
	Real const z1=p2(2);

	if ( std::fabs(x0-x1)<myeps )
		return NAN;

	if (dir==1)
		return y0 + (xval - x0)/(x1 - x0)*(y1 - y0);
	if (dir==2)
		return z0 + (xval - x0)/(x1 - x0)*(z1 - z0);

	exit(1);
	return 0.0;
}

bool doesFaceIntersectElement_3D(RealVectorValue const & p1, RealVectorValue const & p2, RealVectorValue const & p3, RealVectorValue const & hmin, RealVectorValue const & hmax)
{
	RealTensorValue A(	p2(0)-p1(0),p3(0)-p1(0),0.0,
						p2(1)-p1(1),p3(1)-p1(1),0.0,
						p2(2)-p1(2),p3(2)-p1(2),0.0);

	RealVectorValue b=-1.0*p1;

	if ( doesFaceIntersectEdge_3D(A,b,hmin,hmax,0) )
		return true;
	if ( doesFaceIntersectEdge_3D(A,b,hmin,hmax,1) )
		return true;
	if ( doesFaceIntersectEdge_3D(A,b,hmin,hmax,2) )
		return true;

	return false;
}

bool doesFaceIntersectEdge_3D(RealTensorValue const & A, RealVectorValue const & b, RealVectorValue const & hmin, RealVectorValue const & hmax, int const dir)
{
	RealVectorValue const ** ciao;
	ciao = new RealVectorValue const * [2];
	ciao[0]=&hmin;
	ciao[1]=&hmax;
	for (int i=0; i<2;++i)
	{
		for (int j=0; j<2; ++j)
		{
			RealVectorValue bloc=b;
			RealTensorValue Aloc=A;
			if (dir==0)
			{
				bloc(1)=bloc(1)+ciao[i][0](1);
				bloc(2)=bloc(2)+ciao[j][0](2);
			}
			if (dir==1)
			{
				bloc(0)=bloc(0)+ciao[i][0](0);
				bloc(2)=bloc(2)+ciao[j][0](2);
			}
			if (dir==2)
			{
				bloc(0)=bloc(0)+ciao[i][0](0);
				bloc(1)=bloc(1)+ciao[j][0](1);
			}
			bloc(dir)=0.0;
			Aloc(dir,0)=0.0;
			Aloc(dir,1)=0.0;
			Aloc(dir,2)=1.0;
			Real detA=Aloc.det();
			if ( std::fabs(detA)>1e-12 )
			{
				RealVectorValue x;
				Aloc.solve(bloc,x);
				if (0.0-myeps<=x(0) && x(0)<=1.0+myeps && 0.0-myeps<=x(1) && x(1)<=1.0+myeps)
				{
					//std::cout<<x<<std::endl;
					Real otherCoord=A(dir,0)*x(0)+A(dir,1)*x(1)-b(dir);
					if (hmin(dir)-myeps<= otherCoord && otherCoord<= hmax(dir)+myeps )
					{
						delete ciao;
						return true;
					}
				}
			}

		}
	}
	delete ciao;
	return false;
}
