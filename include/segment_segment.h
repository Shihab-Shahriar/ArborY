#pragma once

#include <Kokkos_Core.hpp>
#include "vector.h"

static const float ZERO_TOLERANCE = 1e-06f;


/**
 * @brief Compute the squared distance between two line segments.
 * 
 * This code is adapted from Nvidia Physx
 * at: /geomutils/src/distance/GuDistanceSegmentSegment.cpp
 * 
 * If a segment is represented with two points S0, S1, in the form below:
 * S0 = origin + extent * dir;
 * S1 = origin - extent * dir;
 * 
 * @note If a segment is represented with two points(S0, S1), the overloaded 
 * function below can be used. But this one is recommended for performance.
 * 
 * @tparam PxReal type of floating point
 * 
 * @param origin0 origin of the first segment
 * @param dir0 direction of the first segment
 * @param extent0 extent of the first segment
 * @param origin1 origin of the second segment
 * @param dir1 direction of the second segment
 * @param extent1 extent of the second segment
 * 
 * @return PxReal squared distance between the two segments
*/
template <typename PxReal, typename PxVec3=Vec3<PxReal>>
PxReal SegmentSegmentDistanceSquared(	const PxVec3& origin0, const PxVec3& dir0, PxReal extent0,
											const PxVec3& origin1, const PxVec3& dir1, PxReal extent1)
{
    const PxVec3 kDiff	= origin0 - origin1;
    const PxReal fA01	= -dir0.dot(dir1);
    const PxReal fB0	= kDiff.dot(dir0);
    const PxReal fB1	= -kDiff.dot(dir1);
	const PxReal fC		= kDiff.norm2();
	const PxReal fDet	= Kokkos::fabs(1.0f - fA01*fA01);
    PxReal fS0, fS1, fSqrDist, fExtDet0, fExtDet1, fTmpS0, fTmpS1;

    if (fDet >= ZERO_TOLERANCE)
    {
        // segments are not parallel
        fS0 = fA01*fB1-fB0;
        fS1 = fA01*fB0-fB1;
        fExtDet0 = extent0*fDet;
        fExtDet1 = extent1*fDet;

        if (fS0 >= -fExtDet0)
        {
            if (fS0 <= fExtDet0)
            {
                if (fS1 >= -fExtDet1)
                {
                    if (fS1 <= fExtDet1)  // region 0 (interior)
                    {
                        // minimum at two interior points of 3D lines
                        PxReal fInvDet = 1.0f/fDet;
                        fS0 *= fInvDet;
                        fS1 *= fInvDet;
                        fSqrDist = fS0*(fS0+fA01*fS1+2.0f*fB0) + fS1*(fA01*fS0+fS1+2.0f*fB1)+fC;
                    }
                    else  // region 3 (side)
                    {
                        fS1 = extent1;
                        fTmpS0 = -(fA01*fS1+fB0);
                        if (fTmpS0 < -extent0)
                        {
                            fS0 = -extent0;
                            fSqrDist = fS0*(fS0-2.0f*fTmpS0) + fS1*(fS1+2.0f*fB1)+fC;
                        }
                        else if (fTmpS0 <= extent0)
                        {
                            fS0 = fTmpS0;
                            fSqrDist = -fS0*fS0+fS1*(fS1+2.0f*fB1)+fC;
                        }
                        else
                        {
                            fS0 = extent0;
                            fSqrDist = fS0*(fS0-2.0f*fTmpS0) + fS1*(fS1+2.0f*fB1)+fC;
                        }
                    }
                }
                else  // region 7 (side)
                {
                    fS1 = -extent1;
                    fTmpS0 = -(fA01*fS1+fB0);
                    if (fTmpS0 < -extent0)
                    {
                        fS0 = -extent0;
                        fSqrDist = fS0*(fS0-2.0f*fTmpS0) + fS1*(fS1+2.0f*fB1)+fC;
                    }
                    else if (fTmpS0 <= extent0)
                    {
                        fS0 = fTmpS0;
                        fSqrDist = -fS0*fS0+fS1*(fS1+2.0f*fB1)+fC;
                    }
                    else
                    {
                        fS0 = extent0;
                        fSqrDist = fS0*(fS0-2.0f*fTmpS0) + fS1*(fS1+2.0f*fB1)+fC;
                    }
                }
            }
            else
            {
                if (fS1 >= -fExtDet1)
                {
                    if (fS1 <= fExtDet1)  // region 1 (side)
                    {
                        fS0 = extent0;
                        fTmpS1 = -(fA01*fS0+fB1);
                        if (fTmpS1 < -extent1)
                        {
                            fS1 = -extent1;
                            fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                        }
                        else if (fTmpS1 <= extent1)
                        {
                            fS1 = fTmpS1;
                            fSqrDist = -fS1*fS1+fS0*(fS0+2.0f*fB0)+fC;
                        }
                        else
                        {
                            fS1 = extent1;
                            fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                        }
                    }
                    else  // region 2 (corner)
                    {
                        fS1 = extent1;
                        fTmpS0 = -(fA01*fS1+fB0);
                        if (fTmpS0 < -extent0)
                        {
                            fS0 = -extent0;
                            fSqrDist = fS0*(fS0-2.0f*fTmpS0) + fS1*(fS1+2.0f*fB1)+fC;
                        }
                        else if (fTmpS0 <= extent0)
                        {
                            fS0 = fTmpS0;
                            fSqrDist = -fS0*fS0+fS1*(fS1+2.0f*fB1)+fC;
                        }
                        else
                        {
                            fS0 = extent0;
                            fTmpS1 = -(fA01*fS0+fB1);
                            if (fTmpS1 < -extent1)
                            {
                                fS1 = -extent1;
                                fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                            }
                            else if (fTmpS1 <= extent1)
                            {
                                fS1 = fTmpS1;
                                fSqrDist = -fS1*fS1+fS0*(fS0+2.0f*fB0) + fC;
                            }
                            else
                            {
                                fS1 = extent1;
                                fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                            }
                        }
                    }
                }
                else  // region 8 (corner)
                {
                    fS1 = -extent1;
                    fTmpS0 = -(fA01*fS1+fB0);
                    if (fTmpS0 < -extent0)
                    {
                        fS0 = -extent0;
                        fSqrDist = fS0*(fS0-2.0f*fTmpS0) + fS1*(fS1+2.0f*fB1)+fC;
                    }
                    else if (fTmpS0 <= extent0)
                    {
                        fS0 = fTmpS0;
                        fSqrDist = -fS0*fS0+fS1*(fS1+2.0f*fB1)+fC;
                    }
                    else
                    {
                        fS0 = extent0;
                        fTmpS1 = -(fA01*fS0+fB1);
                        if (fTmpS1 > extent1)
                        {
                            fS1 = extent1;
                            fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                        }
                        else if (fTmpS1 >= -extent1)
                        {
                            fS1 = fTmpS1;
                            fSqrDist = -fS1*fS1+fS0*(fS0+2.0f*fB0) + fC;
                        }
                        else
                        {
                            fS1 = -extent1;
                            fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                        }
                    }
                }
            }
        }
        else 
        {
            if (fS1 >= -fExtDet1)
            {
                if (fS1 <= fExtDet1)  // region 5 (side)
                {
                    fS0 = -extent0;
                    fTmpS1 = -(fA01*fS0+fB1);
                    if (fTmpS1 < -extent1)
                    {
                        fS1 = -extent1;
                        fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                    }
                    else if (fTmpS1 <= extent1)
                    {
                        fS1 = fTmpS1;
                        fSqrDist = -fS1*fS1+fS0*(fS0+2.0f*fB0)+fC;
                    }
                    else
                    {
                        fS1 = extent1;
                        fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                    }
                }
                else  // region 4 (corner)
                {
                    fS1 = extent1;
                    fTmpS0 = -(fA01*fS1+fB0);
                    if (fTmpS0 > extent0)
                    {
                        fS0 = extent0;
                        fSqrDist = fS0*(fS0-2.0f*fTmpS0) + fS1*(fS1+2.0f*fB1)+fC;
                    }
                    else if (fTmpS0 >= -extent0)
                    {
                        fS0 = fTmpS0;
                        fSqrDist = -fS0*fS0+fS1*(fS1+2.0f*fB1)+fC;
                    }
                    else
                    {
                        fS0 = -extent0;
                        fTmpS1 = -(fA01*fS0+fB1);
                        if (fTmpS1 < -extent1)
                        {
                            fS1 = -extent1;
                            fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                        }
                        else if (fTmpS1 <= extent1)
                        {
                            fS1 = fTmpS1;
                            fSqrDist = -fS1*fS1+fS0*(fS0+2.0f*fB0) + fC;
                        }
                        else
                        {
                            fS1 = extent1;
                            fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                        }
                    }
                }
            }
            else   // region 6 (corner)
            {
                fS1 = -extent1;
                fTmpS0 = -(fA01*fS1+fB0);
                if (fTmpS0 > extent0)
                {
                    fS0 = extent0;
                    fSqrDist = fS0*(fS0-2.0f*fTmpS0) + fS1*(fS1+2.0f*fB1)+fC;
                }
                else if (fTmpS0 >= -extent0)
                {
                    fS0 = fTmpS0;
                    fSqrDist = -fS0*fS0+fS1*(fS1+2.0f*fB1)+fC;
                }
                else
                {
                    fS0 = -extent0;
                    fTmpS1 = -(fA01*fS0+fB1);
                    if (fTmpS1 < -extent1)
                    {
                        fS1 = -extent1;
                        fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                    }
                    else if (fTmpS1 <= extent1)
                    {
                        fS1 = fTmpS1;
                        fSqrDist = -fS1*fS1+fS0*(fS0+2.0f*fB0) + fC;
                    }
                    else
                    {
                        fS1 = extent1;
                        fSqrDist = fS1*(fS1-2.0f*fTmpS1) + fS0*(fS0+2.0f*fB0)+fC;
                    }
                }
            }
        }
    }
    else
    {
		// The segments are parallel.
		PxReal fE0pE1 = extent0 + extent1;
		PxReal fSign = (fA01 > 0.0f ? -1.0f : 1.0f);
		PxReal b0Avr = 0.5f*(fB0 - fSign*fB1);
		PxReal fLambda = -b0Avr;
		if(fLambda < -fE0pE1)
		{
			fLambda = -fE0pE1;
		}
		else if(fLambda > fE0pE1)
		{
			fLambda = fE0pE1;
		}

		fS1 = -fSign*fLambda*extent1/fE0pE1;
		fS0 = fLambda + fSign*fS1;
		fSqrDist = fLambda*(fLambda + 2.0f*b0Avr) + fC;
	}

	// account for numerical round-off error
	return Kokkos::fmax(0.0f, fSqrDist);
}


template <typename PxReal, typename PxVec3=Vec3<PxReal>>
PxReal SegmentSegmentDistanceSquared(const PxVec3& origin0, const PxVec3& end0, 
                                const PxVec3& origin1, const PxVec3& end1)
{
    auto extent0 = (end0 - origin0);
    auto extent1 = (end1 - origin1);
    PxVec3 dir0 = extent0;
	const PxVec3 center0 = origin0 + extent0*0.5f;
	PxReal length0 = extent0.norm();	
	const bool b0 = length0 != 0.0f;
	PxReal oneOverLength0 = 0.0f;
	if(b0)
	{
		oneOverLength0 = 1.0f / length0;
		dir0 = dir0 * oneOverLength0;
		length0 *= 0.5f;
	}

	PxVec3 dir1 = extent1;
	const PxVec3 center1 = origin1 + extent1*0.5f;
	PxReal length1 = extent1.norm();
	const bool b1 = length1 != 0.0f;
	PxReal oneOverLength1 = 0.0f;
	if(b1)
	{
		oneOverLength1 = 1.0f / length1;
		dir1 = dir1 * oneOverLength1;
		length1 *= 0.5f;
	}

    // print all parameters before calling
    // std::cout<<"origin0: "<<origin0<<std::endl;
    // std::cout<<"dir0: "<<dir0<<std::endl;
    // std::cout<<"length0: "<<length0<<std::endl;

    // std::cout<<"origin1: "<<origin1<<std::endl;
    // std::cout<<"dir1: "<<dir1<<std::endl;
    // std::cout<<"length1: "<<length1<<std::endl;    

    return SegmentSegmentDistanceSquared<PxReal, PxVec3>(center0, dir0, length0, center1, dir1, length1);
}


/**
 * DONT USE IT.
 * 
 * From geometric Tools. This works. Originally adapted here to investigate
 * CGAL vs Physx discrepancy. Leaving here in case future need arises.
 * 
*/
KOKKOS_INLINE_FUNCTION
double __dist3D_Segment_to_Segment(
    Vec3<double> const& P0, Vec3<double> const& P1,
    Vec3<double> const& Q0, Vec3<double> const& Q1)
{
    double const SMALL_NUM = 0.00000001;
    Vec3<double>   u = P1 - P0;
    Vec3<double>   v = Q1 - Q0;
    Vec3<double>   w = P0 - Q0;
    double    a = u.dot(u);         // always >= 0
    double    b = u.dot(v);
    double    c = v.dot(v);         // always >= 0
    double    d = u.dot(w);
    double    e = v.dot(w);
    double    D = a*c - b*b;        // always >= 0
    double    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
    double    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;         // force using point P0 on segment S1
        sD = 1.0;         // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                 // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (Kokkos::fabs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (Kokkos::fabs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    // get the difference of the two closest points
    Vec3d closest0 = P0 * (1.0 - sc) + P1 * sc;
    Vec3d closest1 = Q0 * (1.0 - tc) + Q1 * tc;
    Vec3d  diff = closest0 - closest1;
    return diff.dot(diff);
}