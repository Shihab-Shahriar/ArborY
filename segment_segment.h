#pragma once

#include <Kokkos_Core.hpp>
#include "vector.h"
#include <iostream>


// Compute the root of h(z) = h0 + slope*z and clamp it to the interval
// [0,1]. It is required that for h1 = h(1), either (h0 < 0 and h1 > 0)
// or (h0 > 0 and h1 < 0).

template <typename T>
KOKKOS_INLINE_FUNCTION
T GetClampedRoot(T const& slope, T const& h0, T const& h1)
{
    // Theoretically, r is in (0,1). However, when the slope is
    // nearly zero, then so are h0 and h1. Significant numerical
    // rounding problems can occur when using floating-point
    // arithmetic. If the rounding causes r to be outside the
    // interval, clamp it. It is possible that r is in (0,1) and has
    // rounding errors, but because h0 and h1 are both nearly zero,
    // the quadratic is nearly constant on (0,1). Any choice of p
    // should not cause undesirable accuracy problems for the final
    // distance computation.
    //
    // NOTE: You can use bisection to recompute the root or even use
    // bisection to compute the root and skip the division. This is
    // generally slower, which might be a problem for high-performance
    // applications.

    T const zero = static_cast<T>(0);
    T r;
    if (h0 < zero)
    {
        T const one = static_cast<T>(1);
        if (h1 > zero)
        {
            r = -h0 / slope;
            if (r > one)
            {
                r = static_cast<T>(0.5);
            }
            // The slope is positive and -h0 is positive, so there is
            // no need to test for a negative value and clamp it.
        }
        else
        {
            r = one;
        }
    }
    else
    {
        r = zero;
    }
    return r;
}

// Compute the intersection of the line dR/ds = 0 with the domain
// [0,1]^2. The direction of the line dR/ds is conjugate to (1,0),
// so the algorithm for minimization is effectively the conjugate
// gradient algorithm for a quadratic function.
template <typename T>
KOKKOS_INLINE_FUNCTION
static void ComputeIntersection(std::array<T, 2> const& sValue,
    std::array<int32_t, 2> const& classify, T const& b, T const& f00,
    T const& f10, std::array<int32_t, 2>& edge,
    std::array<std::array<T, 2>, 2>& end)
{
    // The divisions are theoretically numbers in [0,1]. Numerical
    // rounding errors might cause the result to be outside the
    // interval. When this happens, it must be that both numerator
    // and denominator are nearly zero. The denominator is nearly
    // zero when the segments are nearly perpendicular. The
    // numerator is nearly zero when the P-segment is nearly
    // degenerate (f00 = a is small). The choice of 0.5 should not
    // cause significant accuracy problems.
    //
    // NOTE: You can use bisection to recompute the root or even use
    // bisection to compute the root and skip the division. This is
    // generally slower, which might be a problem for high-performance
    // applications.

    T const zero = static_cast<T>(0);
    T const half = static_cast<T>(0.5);
    T const one = static_cast<T>(1);
    if (classify[0] < 0)
    {
        edge[0] = 0;
        end[0][0] = zero;
        end[0][1] = f00 / b;
        if (end[0][1] < zero || end[0][1] > one)
        {
            end[0][1] = half;
        }

        if (classify[1] == 0)
        {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else  // classify[1] > 0
        {
            edge[1] = 1;
            end[1][0] = one;
            end[1][1] = f10 / b;
            if (end[1][1] < zero || end[1][1] > one)
            {
                end[1][1] = half;
            }
        }
    }
    else if (classify[0] == 0)
    {
        edge[0] = 2;
        end[0][0] = sValue[0];
        end[0][1] = zero;

        if (classify[1] < 0)
        {
            edge[1] = 0;
            end[1][0] = zero;
            end[1][1] = f00 / b;
            if (end[1][1] < zero || end[1][1] > one)
            {
                end[1][1] = half;
            }
        }
        else if (classify[1] == 0)
        {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else
        {
            edge[1] = 1;
            end[1][0] = one;
            end[1][1] = f10 / b;
            if (end[1][1] < zero || end[1][1] > one)
            {
                end[1][1] = half;
            }
        }
    }
    else  // classify[0] > 0
    {
        edge[0] = 1;
        end[0][0] = one;
        end[0][1] = f10 / b;
        if (end[0][1] < zero || end[0][1] > one)
        {
            end[0][1] = half;
        }

        if (classify[1] == 0)
        {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else
        {
            edge[1] = 0;
            end[1][0] = zero;
            end[1][1] = f00 / b;
            if (end[1][1] < zero || end[1][1] > one)
            {
                end[1][1] = half;
            }
        }
    }
}

// Compute the location of the minimum of R on the segment of
// intersection for the line dR/ds = 0 and the domain [0,1]^2.
template <typename T>
KOKKOS_INLINE_FUNCTION
static void ComputeMinimumParameters(std::array<int32_t, 2> const& edge,
    std::array<std::array<T, 2>, 2> const& end, T const& b, T const& c,
    T const& e, T const& g00, T const& g10, T const& g01, T const& g11,
    std::array<T, 2>& parameter)
{
    T const zero = static_cast<T>(0);
    T const one = static_cast<T>(1);
    T delta = end[1][1] - end[0][1];
    T h0 = delta * (-b * end[0][0] + c * end[0][1] - e);
    if (h0 >= zero)
    {
        if (edge[0] == 0)
        {
            parameter[0] = zero;
            parameter[1] = GetClampedRoot(c, g00, g01);
        }
        else if (edge[0] == 1)
        {
            parameter[0] = one;
            parameter[1] = GetClampedRoot(c, g10, g11);
        }
        else
        {
            parameter[0] = end[0][0];
            parameter[1] = end[0][1];
        }
    }
    else
    {
        T h1 = delta * (-b * end[1][0] + c * end[1][1] - e);
        if (h1 <= zero)
        {
            if (edge[1] == 0)
            {
                parameter[0] = zero;
                parameter[1] = GetClampedRoot(c, g00, g01);
            }
            else if (edge[1] == 1)
            {
                parameter[0] = one;
                parameter[1] = GetClampedRoot(c, g10, g11);
            }
            else
            {
                parameter[0] = end[1][0];
                parameter[1] = end[1][1];
            }
        }
        else  // h0 < 0 and h1 > 0
        {
            T z = std::min(std::max(h0 / (h0 - h1), zero), one);
            T omz = one - z;
            parameter[0] = omz * end[0][0] + z * end[1][0];
            parameter[1] = omz * end[0][1] + z * end[1][1];
        }
    }
}

template <typename T>
T SegmentToSegmentDistance(Vec3d const& P0, Vec3d const& P1,
            Vec3d const& Q0, Vec3d const& Q1)
{
    // The code allows degenerate line segments; that is, P0 and P1
    // can be the same point or Q0 and Q1 can be the same point.  The
    // quadratic function for squared distance between the segment is
    //   R(s,t) = a*s^2 - 2*b*s*t + c*t^2 + 2*d*s - 2*e*t + f
    // for (s,t) in [0,1]^2 where
    //   a = Dot(P1-P0,P1-P0), b = Dot(P1-P0,Q1-Q0), c = Dot(Q1-Q0,Q1-Q0),
    //   d = Dot(P1-P0,P0-Q0), e = Dot(Q1-Q0,P0-Q0), f = Dot(P0-Q0,P0-Q0)
    Vec3d P1mP0 = P1 - P0;
    Vec3d Q1mQ0 = Q1 - Q0;
    Vec3d P0mQ0 = P0 - Q0;

    T a = P1mP0.dot(P1mP0);
    T b = P1mP0.dot(Q1mQ0);
    T c = Q1mQ0.dot(Q1mQ0);
    T d = P1mP0.dot(P0mQ0);
    T e = Q1mQ0.dot(P0mQ0);

    // The derivatives dR/ds(i,j) at the four corners of the domain.
    T f00 = d;
    T f10 = f00 + a;
    T f01 = f00 - b;
    T f11 = f10 - b;

    // The derivatives dR/dt(i,j) at the four corners of the domain.
    T g00 = -e;
    T g10 = g00 - b;
    T g01 = g00 + c;
    T g11 = g10 + c;

    T const zero = static_cast<T>(0);
    T const one = static_cast<T>(1);
    T parameter0, parameter1;

    if (a > zero && c > zero)
    {
        // Compute the solutions to dR/ds(s0,0) = 0 and
        // dR/ds(s1,1) = 0.  The location of sI on the s-axis is
        // stored in classifyI (I = 0 or 1).  If sI <= 0, classifyI
        // is -1.  If sI >= 1, classifyI is 1.  If 0 < sI < 1,
        // classifyI is 0.  This information helps determine where to
        // search for the minimum point (s,t).  The fij values are
        // dR/ds(i,j) for i and j in {0,1}.

        // std::array<T, 2> sValue
        // {
        //     GetClampedRoot(a, f00, f10),
        //     GetClampedRoot(a, f01, f11)
        // };

        Vec3d sValue(
            GetClampedRoot(a, f00, f10),
            GetClampedRoot(a, f01, f11),
            0.0
        );


        Vec3<int32_t> classify(
            sValue[0] <= zero ? -1 : (sValue[0] >= one ? +1 : 0),
            sValue[1] <= zero ? -1 : (sValue[1] >= one ? +1 : 0),
            0
        );


        if (classify[0] == -1 && classify[1] == -1)
        {
            // The minimum must occur on s = 0 for 0 <= t <= 1.
            parameter0 = zero;
            parameter1 = GetClampedRoot(c, g00, g01);
        }
        else if (classify[0] == +1 && classify[1] == +1)
        {
            // The minimum must occur on s = 1 for 0 <= t <= 1.
            parameter0 = one;
            parameter1 = GetClampedRoot(c, g10, g11);
        }
        else
        {
            // The line dR/ds = 0 intersects the domain [0,1]^2 in a
            // nondegenerate segment. Compute the endpoints of that
            // segment, end[0] and end[1]. The edge[i] flag tells you
            // on which domain edge end[i] lives: 0 (s=0), 1 (s=1),
            // 2 (t=0), 3 (t=1).
            std::array<int32_t, 2> edge{ 0, 0 };
            std::array<std::array<T, 2>, 2> end{};
            ComputeIntersection(sValue, classify, b, f00, f10, edge, end);

            // The directional derivative of R along the segment of
            // intersection is
            //   H(z) = (end[1][1]-end[1][0]) *
            //          dR/dt((1-z)*end[0] + z*end[1])
            // for z in [0,1]. The formula uses the fact that
            // dR/ds = 0 on the segment. Compute the minimum of
            // H on [0,1].
            ComputeMinimumParameters(edge, end, b, c, e, g00, g10,
                g01, g11, parameter0, parameter1);
        }
    }
    else
    {
        if (a > zero)
        {
            // The Q-segment is degenerate (Q0 and Q1 are the same
            // point) and the quadratic is R(s,0) = a*s^2 + 2*d*s + f
            // and has (half) first derivative F(t) = a*s + d.  The
            // closest P-point is interior to the P-segment when
            // F(0) < 0 and F(1) > 0.
            parameter0 = GetClampedRoot(a, f00, f10);
            parameter1 = zero;
        }
        else if (c > zero)
        {
            // The P-segment is degenerate (P0 and P1 are the same
            // point) and the quadratic is R(0,t) = c*t^2 - 2*e*t + f
            // and has (half) first derivative G(t) = c*t - e.  The
            // closest Q-point is interior to the Q-segment when
            // G(0) < 0 and G(1) > 0.
            parameter0 = zero;
            parameter1 = GetClampedRoot(c, g00, g01);
        }
        else
        {
            // P-segment and Q-segment are degenerate.
            parameter0 = zero;
            parameter1 = zero;
        }
    }


    result.closest[0] =
        (one - parameter0) * P0 + parameter0 * P1;
    result.closest[1] =
        (one - parameter1) * Q0 + parameter1 * Q1;
    Vector<N, T> diff = result.closest[0] - result.closest[1];
    result.sqrDistance = Dot(diff, diff);
    result.distance = std::sqrt(result.sqrDistance);
    return result;
}
