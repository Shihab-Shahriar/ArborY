#pragma once

#include <Kokkos_Core.hpp>
#include "vector.h"


template <typename PxReal, typename PxVec3=Vec3<PxReal>>
void edgeEdgeDist(PxVec3& x, PxVec3& y,				// closest points
				 const PxVec3& p, const PxVec3& a,	// seg 1 origin, vector
				 const PxVec3& q, const PxVec3& b)	// seg 2 origin, vector
{
	const PxVec3 T = q - p;
	const PxReal ADotA = a.dot(a);
	const PxReal BDotB = b.dot(b);
	const PxReal ADotB = a.dot(b);
	const PxReal ADotT = a.dot(T);
	const PxReal BDotT = b.dot(T);

	// t parameterizes ray (p, a)
	// u parameterizes ray (q, b)

	// Compute t for the closest point on ray (p, a) to ray (q, b)
	const PxReal Denom = ADotA*BDotB - ADotB*ADotB;

	PxReal t;	// We will clamp result so t is on the segment (p, a)
	if(Denom!=0.0f)	
		t = Kokkos::clamp((ADotT*BDotB - BDotT*ADotB) / Denom, 0.0, 1.0);
	else
		t = 0.0f;

	// find u for point on ray (q, b) closest to point at t
	PxReal u;
	if(BDotB!=0.0f)
	{
		u = (t*ADotB - BDotT) / BDotB;

		// if u is on segment (q, b), t and u correspond to closest points, otherwise, clamp u, recompute and clamp t
		if(u<0.0f)
		{
			u = 0.0f;
			if(ADotA!=0.0f)
				t = Kokkos::clamp(ADotT / ADotA, 0.0, 1.0);
			else
				t = 0.0f;
		}
		else if(u > 1.0f)
		{
			u = 1.0f;
			if(ADotA!=0.0f)
				t = Kokkos::clamp((ADotB + ADotT) / ADotA, 0.0, 1.0);
			else
				t = 0.0f;
		}
	}
	else
	{
		u = 0.0f;
		if(ADotA!=0.0f)
			t = Kokkos::clamp(ADotT / ADotA, 0.0, 1.0);
		else
			t = 0.0f;
	}

	x = p + a * t;
	y = q + b * u;
}


/**
 * The code has been adapted from Nvidia Physx library. Link:
 * /source/geomutils/src/distance/GuDistanceTriangleTriangle.cpp
 * 
 * The code doesn't have any comment. To understand and annotate it, 
 * this SO answer was helpful: https://stackoverflow.com/a/53622315/4553309
 * 
 * @tparam PxReal 
 * @tparam PxVec3
 * 
 * @param p The vertices of first triangle
 * @param q The vertices of second triangle
 * 
 * @return PxReal The squared distance between the two triangles
 * 
*/
template <typename PxReal, typename PxVec3=Vec3<PxReal>>
float TriangleTriangleDistanceSquared(const PxVec3 p[3], const PxVec3 q[3])
{
	PxVec3 Sv[3];
    Sv[0] = p[1] - p[0];
    Sv[1] = p[2] - p[1];
    Sv[2] = p[0] - p[2];

	PxVec3 Tv[3];
    Tv[0] = q[1] - q[0];
    Tv[1] = q[2] - q[1];
    Tv[2] = q[0] - q[2];

	PxVec3 minP, minQ;
    PxVec3 cp, cq; // the closest points from the two triangles
	bool shown_disjoint = false; // are 2 trianges disjoint?

	float mindd = Kokkos::Experimental::finite_max_v<float>;;

	// find the minimum distance between each pair of edges
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			edgeEdgeDist<PxReal, PxVec3>(cp, cq, p[i], Sv[i], q[j], Tv[j]);
			const PxVec3 V = cq - cp;
			const float dd = V.dot(V);

			// new min distance edge pair
			if(dd<=mindd)
			{
				minP = cp;
				minQ = cq;
				mindd = dd;

				//--------------------------
				// Checks if off-edge vertices are on the opposite side of 
				// plane stab parallel to the vector V (connecting the closest points)
				// If so, we can be certain we have found nearest points, and return.
				int id = i+2;
				if(id>=3)
					id-=3;
				PxVec3 Z = p[id] - cp;
				float a = Z.dot(V);
				id = j+2;
				if(id>=3)
					id-=3;
				Z = q[id] - cq;
				float b = Z.dot(V);

				if((a<=0.0f) && (b>=0.0f)){
					return V.dot(V);
				}
				//----END--------------

				if(a<=0.0f)	a = 0.0f;
				else if(b>0.0f)	b = 0.0f;

				if((mindd - a + b) > 0.0f)
					shown_disjoint = true;
			}
		}
	}

	// std::cout<<"Didn't return early"<<std::endl;
	// std::cout<<cq<<std::endl;

	/**
	 * If we're here, one of 4 edge cases to consider now:
	 * 1. One of the closest points is a vertex of one triangle and the other 
	 * closest point is on the face of the other triangle 
	 * 2. The triangles intersect 
	 * 3. An edge of one triangle is parallel to the face of the other triangle 
	 * 4. One or both triangles are degenerate
	 * 
	 * The big chuck of code below is mostly for edge case 1. For 2 and 4, the
	 * previously found closest points are good enough.
	*/
	PxVec3 Sn = Sv[0].cross(Sv[1]);
	float Snl = Sn.dot(Sn);

	// triangle is not degenerate
	if(Snl>1e-15f)
	{
		const PxVec3 Tp((p[0] - q[0]).dot(Sn),
						(p[0] - q[1]).dot(Sn),
						(p[0] - q[2]).dot(Sn));

		// If all vertices of triangle q are on the same side of plane of 
		// triangle p, the signs should be same. Fidn the closest vertex
		// on variable `index`.
		int index = -1;
		if((Tp[0]>0.0f) && (Tp[1]>0.0f) && (Tp[2]>0.0f))
		{
			if(Tp[0]<Tp[1])		index = 0; else index = 1;
			if(Tp[2]<Tp[index])	index = 2;
		}
		else if((Tp[0]<0.0f) && (Tp[1]<0.0f) && (Tp[2]<0.0f))
		{
			if(Tp[0]>Tp[1])		index = 0; else index = 1;
			if(Tp[2]>Tp[index])	index = 2;
		}

		if(index >= 0) // All on same side, proceed with edge case 1
		{
			shown_disjoint = true;

			const PxVec3& qIndex = q[index];

			// The 3 'if's belows check if the projection of the closest 
			// vertex is inside the triangle p.
			PxVec3 V = qIndex - p[0];
			PxVec3 Z = Sn.cross(Sv[0]);
			if(V.dot(Z)>0.0f)
			{
				V = qIndex - p[1];
				Z = Sn.cross(Sv[1]);
				if(V.dot(Z)>0.0f)
				{
					V = qIndex - p[2];
					Z = Sn.cross(Sv[2]);
					if(V.dot(Z)>0.0f)
					{
						cp = qIndex + Sn * Tp[index] *(1/Snl);
						cq = qIndex;
						return (cp - cq).norm2();
					}
				}
			}
		}
	}

	// repeat above from triangle q's point of view

	PxVec3 Tn = Tv[0].cross(Tv[1]);
	float Tnl = Tn.dot(Tn);
  
	if(Tnl>1e-15f)
	{
		const PxVec3 Sp((q[0] - p[0]).dot(Tn),
						(q[0] - p[1]).dot(Tn),
						(q[0] - p[2]).dot(Tn));

		int index = -1;
		if((Sp[0]>0.0f) && (Sp[1]>0.0f) && (Sp[2]>0.0f))
		{
			if(Sp[0]<Sp[1])		index = 0; else index = 1;
			if(Sp[2]<Sp[index])	index = 2;
		}
		else if((Sp[0]<0.0f) && (Sp[1]<0.0f) && (Sp[2]<0.0f))
		{
			if(Sp[0]>Sp[1])		index = 0; else index = 1;
			if(Sp[2]>Sp[index])	index = 2;
		}

		if(index >= 0)
		{ 
			shown_disjoint = true;

			const PxVec3& pIndex = p[index];

			PxVec3 V = pIndex - q[0];
			PxVec3 Z = Tn.cross(Tv[0]);
			if(V.dot(Z)>0.0f)
			{
				V = pIndex - q[1];
				Z = Tn.cross(Tv[1]);
				if(V.dot(Z)>0.0f)
				{
					V = pIndex - q[2];
					Z = Tn.cross(Tv[2]);
					if(V.dot(Z)>0.0f)
					{
						cp = pIndex;
						cq = pIndex + Tn * Sp[index]* (1/Tnl);
						return (cp - cq).norm2();
					}
				}
			}
		}
	}

	if(shown_disjoint) return mindd;
	else return 0.0f;
}
