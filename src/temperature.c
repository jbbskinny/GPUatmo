/*
 * File:   temperature.c
 * Author: Justin Blythe
 *
 * Created on August 29, 2012, 7:35 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "radtransfer.h"
#include "temperature.h"
#include "extinction.h"
#include "grids.h"
#include "numerical.h"
#include "utilities.h"

/*****************************************************************************
 *************************** Temperature Profile *****************************
 *****************************************************************************/

void rosseland_temp_profile(const int D, const int K, const int maxiteration, 
                            const double g, const double teff, double *E1, 
                            double *h, double *rho, double *T6, double *tauTh)
/**
 * \brief
 *
 * This function sets up the Rosseland (initial) temperature profile.
 *
 * \param[in]  D            the total number of depth grid points
 * \param[in]  K            the total number of energy grid points
 * \param[in]  maxiteration the maximum number of iterations for convergent 
 *                          temperature profile
 * \param[in]  g            the gravitational acceleration 
 *                          (units: 10^14 cm/s)
 * \param[in]  teff         the effective temperature of the neutron star 
 *                          (units: 10^6 K)
 * \param[in]  E1           the energy grid (units: keV)
 * \param[in]  h            the integration weights of the energy grid
 *                          (units: keV)
 * \param[out] rho          the density profile (units: g cm^-3)
 * \param[out] T6           the temperature profile (units: 10^6 K)
 * \param[in]  tauTh        the Thompson depth grid (units: none)
 *
 * \warning The if statement towards the end of this function makes sure that 
 *          the temperature profile converged within the alotted number of 
 *          iterations.  If it does not, it exits the program and prints an 
 *          error message.
 */                            
{
    update_rossld_profile(D, teff, T6, tauTh);
    update_density_profile(D, g, rho, T6, tauTh);

    double *oldT6 = malloc_double_1d(D);
    double *chiR  = malloc_double_1d(D); 
    double *tauR  = malloc_double_1d(D);

    int i;

    for(i = 0; i < maxiteration; i++)
    {
        calc_rossld_extinc(D, K, chiR, E1, h, rho, T6);
        calc_rossld_depth(D, chiR, tauR, tauTh);
        update_rossld_profile(D, teff, T6, tauR);
        update_density_profile(D, g, rho, T6, tauTh);

        if(fabs(residual(D, 1e-4, 0.0, oldT6, T6)) < 1.0) break;

        copy_array(D, oldT6, T6);
    }
    
    if(i == maxiteration)
    {
        printf("ERROR in %s at line %d\n", __FILE__, __LINE__);
        exit(-1);
    }

    free_1d(oldT6); 
    free_1d(chiR); 
    free_1d(tauR);
}

void update_rossld_profile(const int D, const double teff, double *T6, 
                           double *tau)
/**
 * \brief
 *
 * This function updates the temperature profile based on the current values of
 * the Rosseland depth. The temperature profile is initialized using Thompson 
 * depth.
 *
 * \param[in]  D    the total number of depth grid points
 * \param[in]  teff the effective temperature of the neutron star
 *                  (units: 10^6 K)
 * \param[out] T6   the temperature profile (units: 10^6 K)
 * \param[in]  tau  the depth grid, either Thompson (initially) or Rosseland
 */
{
    for(int d = 0; d < D; d++) T6[d] = teff * pow((0.75 * tau[d] + 0.5), 0.25);
}

void calc_rossld_depth(const int D, double *chiR, double *tauR, double *tauTh)
/**
 * \brief
 * 
 * This function is used to calculate the Rosseland depth.
 *
 * \param[in]  D     the total number of depth grid points
 * \param[in]  chiR  the Rosseland extinction coefficient
 * \param[out] tauR  the Rosseland depth grid
 * \param[in]  tauTh the Thompson depth grid
 */
{
    double delta = 1.15130 * (log10(tauTh[1]) - log10(tauTh[0]));
    
    tauR[0]      = tauTh[0]; 

    for(int d = 1; d < D; d++) tauR[d] = tauR[d-1] + (tauTh[d] * chiR[d] + 
                                         tauTh[d-1] * chiR[d-1]) * delta;
}

void calc_rossld_extinc(const int D, const int K, double *chiR, double *E1, 
                        double *h, double *rho, double *T6)
/**
 * \brief
 *
 * This function is used to calculate the Rosseland extinction coefficient at 
 * at each depth.
 *
 * \param[in]  D    the total number of depth grid points
 * \param[in]  K    the total number of energy grid points
 * \param[out] chiR the Rosseland extinction coefficient
 * \param[in]  E1   the energy grid (units: keV)
 * \param[in]  h    the integration weights of the energy grid
 *                  (units: keV)
 * \param[in]  rho  the density profile (units: g cm^-3)
 * \param[in]  T6   the temperature profile (units: 10^6 K)
 */
{
    double *f = malloc_double_1d(K);
    
    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++)
        {
            double chi   = chi_ff(E1[k], rho[d], T6[d]);
            double boltz = exp(11.60445 * E1[k] / T6[d]);

            f[k]         = 8.101324e3 * pow(E1[k], 4.0) * boltz /
                           (pow(boltz - 1.0, 2.0) * pow(T6[d], 5.0) * chi);
        }

        chiR[d] = 1 / (quadrature(K, f, h));
    }

    free(f);
}
