/*
 * File:   numerical.c
 * Author: Justin Blythe
 *
 * Created on September 10, 2012, 1:18 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*****************************************************************************
 ************************** Numerical Routines *******************************
 *****************************************************************************/

double clenshaw(const int N, double x, double *c)
/**
 * \brief
 *
 * This function uses Clenshaw's recurrance to get an approximate value to a
 * function, y = f(x).  It takes in Chebyshev coefficients to a particular 
 * function and an x value to approximate the function f(x) at x.
 *
 * \param[in] N the maximum degree of the polynomial
 * \param[in] x the point to evaluate the function
 * \param[in] c the chebyshev polynomial coefficients of a given function
 *
 * \return the approximate value of the function at x
 */
{
    double dummy, d = 0.0, dd = 0.0;
    double x2 = 2 * x;

    for(int i = N - 1; i >= 1; i--)
        {
            dummy = d;
            d     = x2 * d - dd + c[i];
            dd    = dummy;
        }

    double y = (x * d) - dd + c[0];

    return(y);
}

void setup_gauss_legendre(const int N, const double init, const double final,
                          double *abscis, double *weight)
/**
 * \brief
 *
 * This function calculates the abscissas and weights of the Gauss-Legendre 
 * N-point quadrature formula.  It takes in the limits of integration and N, to
 * calculate the abscis and weights.
 *
 * \param[in]  N      the number of points within the quadrature
 * \param[in]  init   the lower limit of integration
 * \param[in]  final  the upper limit of integration
 * \param[out] abscis the point where to evalulate a function
 * \param[out] weight the integration weight corresponding to the abscis
 */
{
    int    m  = (N + 1) / 2;
    double xm = 0.5 * (final + init);
    double xl = 0.5 * (final - init);

    for(int i = 1; i <= m; i++)
    {
        double pp, z1;
        double z = cos(3.14159 * (i - 0.25) / (N + 0.5));

        do 
        {
            double p1 = 1.0; 
            double p2 = 0.0;
            double p3;

            for(int j = 1; j <= N; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1) * z * p2 - (j - 1.0) * p3) / j;
            }

            pp = N * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z  = z1 - (p1 / pp);

        } while (fabs(z - z1) > 3.0e-11);
        
        abscis[i-1] = xm - (xl * z);
        abscis[N-i] = xm + (xl * z);
        weight[i-1] = (2.0 * xl) / ((1.0 - (z * z)) * (pp * pp));
        weight[N-i] = weight[i-1];
    }
}

double quadrature(const int N, double *f, double *weight)
/**
 * \brief
 *
 * This function evaluates an integral. It takes in the number of points, the 
 * integrand at a given point, and the weight at that point to return the value
 * of the integral.
 *
 * \param[in] N      the number of points
 * \param[in] f      the integrand, length N
 * \param[in] weight the integration weight
 *
 * \return the value of the integral
 */
{
    double out = 0.0;

    for(int i = 0; i < N; i++) out += weight[i] * f[i];

    return(out);
}

