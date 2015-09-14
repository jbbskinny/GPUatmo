/*
 * File: grids.c
 * Author: Justin Blythe
 *
 * Created on August 29, 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grids.h"
#include "numerical.h"
#include "radtransfer.h"
#include "extinction.h"
#include "utilities.h"

/*****************************************************************************
 ********************************** Grids ************************************
 *****************************************************************************/

void setup_tauTh_grid(const int decades, const int ppd, double init, 
                      double *tauTh)
/**
 * \brief
 *
 * This function initializes the Thompson depth grid, given the number of 
 * decades, the points per decade, and the initial value at the top of the 
 * atmosphere. The grid is logarithmically spaced so the sampling within each 
 * order of magnitude is the same.
 *
 * \param[in]  decades the number of orders of magnitude to define tauTh over
 * \param[in]  ppd     the number of grid points within each decade
 * \param[in]  init    the initial value of tauTh, corresponding to the top of 
 *                     the neutron star atmosphere
 * \param[out] tauTh   the Thompson depth grid, logarithmically spaced
 */
{
    int D        = ppd * decades;
    double delta = (double)decades / (D - 1);
    
    tauTh[0] = init;

    for(int d = 1; d < D; d++) tauTh[d] = tauTh[d-1] * pow(10.0, delta);
}

void setup_E1_grid(const int decades, const int ppd, double init, double *E1,
                   double *h)
/**
 * \brief
 *
 * This function initializes the energy grid and it's integration weights, 
 * given the number of decades, the points per decade, and the lowest energy of
 * a photon that one would find in an enviroment with these temperatures.  The 
 * grid is logarithmically spaced so sampling within each order of magnitude is
 * consistent.
 *
 * \param[in]  decades the number of orders of magnitude to define E1 over
 * \param[in]  ppd     the number of grid points within each decade
 * \param[in]  init    the initial value of E1
 * \param[out] E1      the energy grid, logarithmically spaced
 * \param[out] h       the integration weights of the energy grid
 */
{
    int K        = ppd * decades;
    double delta = (double)decades / (K - 1);
    
    E1[0] = init;
    h [0] = 1.15129 * E1[0] * delta;

    for(int k = 1; k < K; k++)
    {
        E1[k] = E1[k-1] * pow(10.0, delta);
        h[k]  = 2.30259 * E1[k]   * delta;

        if(k == K - 1) h[k] = 0.5 * h[k];
    }
}

void setup_mu_grid(const int M, double *mu, double *w)
/**
 * \brief
 * 
 * This function initiallizes the angular grid, more specifically the cosine of
 * the angle. Using Gauss-Legendre quadrature to define the angular grid and 
 * it's integration weights. Since we are taking the cosine of the angle the 
 * grid is defined between 0 and 1. This algorithm is used so that certain, 
 * propbable angles have higher sampling around that angle.
 *
 * \param[in]  M  the total number of angular points
 * \param[out] mu the angular dependent grid (units: None)
 * \param[out] w  the integration weights of the angular dependent grid
 *                (units: None)
 */
{
    setup_gauss_legendre(M, 0.0, 1.0, mu, w);
}

void update_density_profile(const int D, const double g, double *rho,
                            double *T6, double *tauTh)
/**
 * \brief
 *
 * This function updates the density profile based on the current values of the
 * temperature profile.
 *
 * \param[in]  D     the total number of depth grid points
 * \param[in]  g     the gravitational acceleration (units: 10^14 dynes)
 * \param[out] rho   the density profile (units: g cm^-3)
 * \param[in]  T6    the temperature profile (units: 10^6 K)
 * \param[in]  tauTh the Thomson depth grid
 */                            
{
    for(int d = 0; d < D; d++) rho[d] = 1.52297 * g * tauTh[d] / T6[d];
}
