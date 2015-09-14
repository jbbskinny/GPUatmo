/**
 * File: spectral.c
 * Author: Justin Blythe
 *
 * Created on March 11, 2013, 8:30 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "extinction.h" 
#include "numerical.h"
#include "radtransfer.h"
#include "utilities.h"

/*****************************************************************************
 ************************** Spectral Quantities ******************************
 *****************************************************************************/

void update_jnu_grid(const int D, const int K, const int M, double *w,
                     double **iminus, double **iplus, double **jnu)
/**
 * \brief
 *
 * This function calculates the spectral density per unit photon energy in the
 * neutron star atmosphere (NSA). It averages the iplus and iminus beams over
 * all angles to find the amount of intensity at each depth and photon energy.
 *
 * \param[in]  D      the total number of depth grid points
 * \param[in]  K      the total number of energy grid points
 * \param[in]  M      the total number of angular grid points
 * \param[in]  w      the integration weights of the angular dependent grid
 * \param[in]  iminus the specific intensity in the minus z direction
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1 sr^-1)
 * \param[in]  iplus  the specific intensity in the plus z direction
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1 sr^-1)
 * \param[out] jnu    the spectral density per unit photon energy
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1)
 */
{
    double **f = malloc_double_2d(K, M);

    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++)
        {
            for(int m = 0; m < M; m++) f[k][m] = 0.5 * (iplus[d][K*m+k] + 
                                                 iminus[d][K*m+k]);

            jnu[d][k] = quadrature(M, f[k], w);
        }
    }
    
    free_2d(K, f);
}

void update_j_profile(const int D, const int K, double *h, double *j,
                      double **jnu)
/**
 * \brief
 *
 * This function calculates the spectral density in the NSA at all depths.
 *
 * \param[in]  D   the total number of depth grid points
 * \param[in]  K   the total number of energy grid points
 * \param[in]  h   the integration weights of the energy grid
 * \param[out] j   the spectral density of the NSA
 *                 (units: 10^22 erg cm^-2 s^-1)
 * \param[in]  jnu the spectral density per unit photon energy in the NSA
 *                 (units: 10^22 erg cm^-2 s^-1 keV-1)
 */
{
    for(int d = 0; d < D; d++) j[d] = quadrature(K, jnu[d], h);
}

void update_chiJ_profile(const int D, const int K, double *chiJ, double *h, 
                         double *j, double **chi, double **jnu)
/**
 * \brief
 *
 * This function calculates average extinction coefficient, wrt spectral
 * density.
 *
 * \param[in]  D    the total number of depth grid points
 * \param[in]  K    the total number of energy grid points
 * \param[out] chiJ the average extinction wrt spectral density
 *                  (units: None)
 * \param[in]  h    the integration weights of the energy grid
 * \param[in]  j    the spectral density of the NSA
 *                  (units: 10^22 erg cm^-2 s^-1)
 * \param[in]  chi  the extinction coefficient (units: None)
 * \param[in]  jnu  the spectral density per unit photon energy in the NSA
 *                  (units: 10^22 erg cm^-2 s^-1 keV-1)
 */
{
    double **f = malloc_double_2d(D, K);

    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++) f[d][k] = chi[d][k] * jnu[d][k];
        
        chiJ[d] = quadrature(K, f[d], h) / j[d];
    }

    free_2d(D, f);
}

void update_fnu_grid(const int D, const int K, const int M, double *mu, 
                     double *w, double **fnu, double **iminus, double **iplus)
/**
 * \brief
 *
 * This function calculates the spectral flux per unit photon energy in the
 * neutron star atmosphere (NSA). It averages the iplus and iminus beams over
 * all angles to find the amount of intensity moving in the +z direction at 
 * each depth and photon energy.
 *
 * \param[in]  D      the total number of depth grid points
 * \param[in]  K      the total number of energy grid points
 * \param[in]  M      the total number of angular grid points
 * \param[in]  mu     the angular dependent grid (units: None)
 * \param[in]  w      the integration weights of the angular dependent grid
 * \param[out] fnu    the spectral flux per unit photon energy
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1)
 * \param[in]  iminus the specific intensity in the -z direction
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1 sr^-1)
 * \param[in]  iplus  the specific intensity in the +z direction
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1 sr^-1)
 */
{
    double **f = malloc_double_2d(K, M);

    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++)
        {
            for(int m = 0; m < M; m++) f[k][m] = 0.5 * mu[m] * 
                                                 (iplus[d][K*m+k]
                                                 - iminus[d][K*m+k]);

            fnu[d][k] = quadrature(M, f[k], w);
        }
    }

    free_2d(K, f);
}

void update_fbar_profile(const int D, const int K, double *fbar, double *h, 
                         double **fnu)
/**
 * \brief
 *
 * This function calculates the spectral flux in the NSA at all depths.
 *
 * \param[in]  D    the total number of depth grid points
 * \param[in]  K    the total number of energy grid points
 * \param[out] fbar the spectral flux of the NSA
 *                  (units: 10^22 erg cm^-2 s^-1)
 * \param[in]  h    the integration weights of the energy grid
 * \param[in]  jnu  the spectral flux per unit photon energy in the NSA
 *                  (units: 10^22 erg cm^-2 s^-1 keV-1)
 */
{
    for(int d = 0; d < D; d++) fbar[d] = quadrature(K, fnu[d], h);
}

void update_chiF_profile(const int D, const int K, double *chiF, double *fbar,
                         double *h, double **chi, double **fnu)
/**
 * \brief
 *
 * This function calculates average extinction coefficient, wrt spectral
 * flux.
 *
 * \param[in]  D    the total number of depth grid points
 * \param[in]  K    the total number of energy grid points
 * \param[out] chiF the average extinction wrt spectral flux
 *                  (units: None)
 * \param[in]  fbar the spectral flux of the NSA
 *                  (units: 10^22 erg cm^-2 s^-1)
 * \param[in]  h    the integration weights of the energy grid
 * \param[in]  chi  the extinction coefficient (units: None)
 * \param[in]  fnu  the spectral flux per unit photon energy in the NSA
 *                  (units: 10^22 erg cm^-2 s^-1 keV-1)
 */
{
    double **f = malloc_double_2d(D, K);

    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++) f[d][k] = chi[d][k] * fnu[d][k];

        chiF[d] = quadrature(K, f[d], h) / fbar[d];
    }

    free_2d(D, f);
}

void update_b_profile(const int D, double *b, double *T6)
/**
 * \brief
 *
 * This function calculates the photon energy integrated Planck function.
 * Given by an analytic formula.
 *
 * \param[in]  D  the total number of depth grid points
 * \param[out] b  the photon energy intengrated Planck function
 *                (units: 10^22 erg cm^2 s^-1)
 * \param[in]  T6 the temperature profile (units: 10^6 K)
 */
{
    for(int d = 0; d < D; d++) b[d] = 1.80498e-3 * pow(T6[d], 4.0);
}

void update_chiP_profile(const int D, const int K, double *b, double *chiP,
                         double *h, double **chi, double **source)
/**
 * \brief
 *
 * This function calculates the average extinction coefficient, wrt the
 * Planck function.
 *
 * \param[in]  D      the total number of depth grid points
 * \param[in]  K      the total number of energy grid points
 * \param[in]  b      the photon energy intengrated Planck function
 *                    (units: 10^22 erg cm^2 s^-1)
 * \param[out] chiP   the average extinction wrt Planck function
 *                    (units: None)
 * \param[in]  h      the integration weights of the energy grid
 * \param[in]  chi    the extinction coefficient (units: None)
 * \param[in]  source the source function
 *                    (units: 10^22 erg cm^-2 s^-1 keV^-1 sr^-1)
 */
{
    double **f = malloc_double_2d(D, K);

    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++) f[d][k] = chi[d][k] * source[d][k];
        
        chiP[d] = quadrature(K, f[d], h) / b[d];
    }

    free_2d(D, f);
}
