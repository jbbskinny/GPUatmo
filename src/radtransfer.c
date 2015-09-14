/**
 * File: radtransfer.c
 * Author: Justin Blythe
 *
 * Created on October 4, 2012, 11:27 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical.h"
#include "extinction.h" 
#include "radtransfer.h"
#include "utilities.h"

/*****************************************************************************
 ****************** Radiative Transfer Equation Solver ***********************
 *****************************************************************************/

void solve_rte(const int D, const int I, const int K, const int M, double teff,
               double *E1, double *h, double *mu, double *T6, double **chi,
               double **dtau, double **dtaum, double **dtaup, double **iminus,
               double **iplus, double **source)
/**
 * \brief 
 *
 * This function solves the radiative transfer equation, calculating the 
 * specific intensity as a beam of light from core travels through the neutron
 * star atmosphere. The core acts as a perfect blackbody, thus the inner
 * boundary condition is Planck's spectrum (source[D-1]) and the outer 
 * boundary is set to zero.
 *
 * \param[in]  D      the total number of depth grid points
 * \param[in]  I      the total number of energy-angular grid points combined,
 *                    I = K * M
 * \param[in]  K      the total number of energy grid points
 * \param[in]  M      the total number of angular grid points
 * \param[in]  teff   the effective temperature (units: 10^6 K)
 * \param[in]  E1     the energy grid (units: keV)
 * \param[in]  h      the integration weights of the energy grid (units: keV)
 * \param[in]  mu     the angular dependent grid (units: None)
 * \param[in]  T6     the temperature profile (units: 10^6 K)
 * \param[in]  chi    the extinction grid (units: None)
 * \param[in]  dtau   the sum of dtaup and dtaum divided by 2, 2D with D-row
 *                    and I-column
 * \param[in]  dtaum  the optical depth difference defined for the minus 
 *                    directed beam, 2D with D-row and I-column
 * \param[in]  dtaup  the optical depth difference defined for the plus 
                      directed beam, 2D with D-row and I-column
 * \param[out] iminus the specific intensity in the minus z direction
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1 sr^-1)
 * \param[out] iplus  the specific intensity in the plus z direction
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1 sr^-1)
 * \param[in]  source the source function
 *                    (units: 10^22 erg cm^-2 s^-1 keV-1 sr^-1)
 */                     
{
    for(int d = D - 1; d >= 0; d--)
    {
        for(int i = 0; i < I; i++)
        {
            if(d == 0)
            {
                double alphap = 1.0 - exp(-dtaup[0][i]);
                double betap  = alphap - dtaup[0][i] * 
                                exp(-dtaup[0][i]);
                double dsrcp  = (source[1][i] - source[0][i]) / 
                                dtaup[0][i];

                iplus[0][i] = iplus[1][i] * exp(-dtaup[0][i]) + 
                                  alphap * source[0][i] + betap * dsrcp;
                continue;
            }

            if(d == D-1)
            {
                iplus[d][i] = source[D-1][i];
                continue;
            }

            double fp     = dtaum[d][i] / (2.0 * dtaup[d][i] * 
                            dtau[d][i]);
            double fm     = dtaup[d][i] / (2.0 * dtaum[d][i] * 
                            dtau[d][i]);
            double gp     = 1.0 / (dtaup[d][i] * dtau[d][i]);
            double gm     = 1.0 / (dtaum[d][i] * dtau[d][i]);

            double alphap = 1.0 - exp(-dtaup[d][i]);
            double betap  = -dtaup[d][i] * exp(-dtaup[d][i]) + alphap;
            double gammap = -pow(dtaup[d][i], 2.0) * exp(-dtaup[d][i])
                            + 2.0 * betap;

            double ap     = -betap * fm + 0.5 * gammap * gm;
            double bp     = alphap - betap * (fp - fm) - 0.5 * gammap * 
                            (gp + gm);
            double cp     = betap * fp + 0.5 * gammap * gp;

            iplus[d][i]   = iplus[d+1][i] * exp(-dtaup[d][i]) + 
                            ap * source[d-1][i] + bp * source[d][i] + 
                            cp * source[d+1][i];
        }
    }
    
    double *dBdtau = malloc_double_1d(I);

    calc_dBdtau_grid(D, K, M, teff, dBdtau, E1, h, mu, T6, chi);

    for(int d = 0; d < D; d++)
    {
        for(int i = 0; i < I; i++)
        {
            if(d == 0)
            {
                iminus[d][i] = 1e-15;
                continue;
            }

            if(d == D-1)
            {
                iminus[d][i] = iplus[d][i] - dBdtau[i];

                continue;
            }

            double fp     = dtaum[d][i] / (2.0 * dtaup[d][i] * dtau[d][i]);
            double fm     = dtaup[d][i] / (2.0 * dtaum[d][i] * dtau[d][i]);
            double gp     = 1.0 / (dtaup[d][i] * dtau[d][i]);
            double gm     = 1.0 / (dtaum[d][i] * dtau[d][i]);

            double alpham = 1.0 - exp(-dtaum[d][i]);
            double betam  = -dtaum[d][i] * exp(-dtaum[d][i]) + alpham;
            double gammam = -pow(dtaum[d][i], 2.0) * exp(-dtaum[d][i]) + 
                            2.0 * betam;

            double am     = betam * fm + 0.5 * gammam * gm;
            double bm     = alpham + betam * (fp - fm) - 0.5 * gammam * 
                            (gp + gm);
            double cm     = -betam * fp + 0.5 * gammam * gp;

            iminus[d][i] = iminus[d-1][i] * exp(-dtaum[d][i]) + am * 
                           source[d-1][i] + bm * source[d][i] + 
                           cm * source[d+1][i];
        }
    }

    free_1d(dBdtau);
}

void update_dtau_grids(const int D, const int K, const int M, double *mu, 
                       double *tauTh, double **chi, double **dtau, 
                       double **dtaum, double **dtaup)
/**
 * \brief
 *
 * This function calculates the change in the actual optical depth between 
 * subsequent grid points. These arrays are then used by solve_rte() when 
 * calculating the specific intensity profiles.
 *
 * \param[in]  D     the total number of depth grid points
 * \param[in]  K     the total number of energy  grid points
 * \param[in]  M     the total number of angular dependent grid points
 * \param[in]  mu    the angle grid, length M
 * \param[in]  tauTh the Thompson depth grid, length D
 * \param[in]  chi   the extinction coefficient (units: None)
 * \param[out] dtau  the sum of dtaup and dtaum divided by 2, 2D with D-row
 *                   and I-column
 * \param[out] dtaum the optical depth difference defined for the minus 
 *                   directed beam, 2D with D-row and I-column
 * \param[out] dtaup the optical depth difference defined for the plus directed
 *                   beam, 2D with D-row and I-column
 */
{
    double delta  = 1.15129 * (log10(tauTh[1]) - log10(tauTh[0]));

    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++)
        {
            for(int m = 0; m < M; m++)
            {
                if(d == 0) 
                {
                    dtaup[d][K*m+k] = (tauTh[d] * chi[d][k] + 
                                      tauTh[d+1] * chi[d+1][k]) * 
                                      delta / mu[m];
                    continue;
                }

                dtaum[d][K*m+k] = (tauTh[d] * chi[d][k] + 
                                  tauTh[d-1] * chi[d-1][k]) * 
                                  delta / mu[m];

                if(d == D - 1) continue;

                dtaup[d][K*m+k] = (tauTh[d] * chi[d][k] + 
                               tauTh[d+1] * chi[d+1][k]) * 
                               delta / mu[m];

                //if(dtaum[d][K*m+k] < 1e-11) dtaum[d][K*m+k] = 1e-11;
                //if(dtaup[d][K*m+k] < 1e-11) dtaup[d][K*m+k] = 1e-11;


                dtau[d][K*m+k] = 0.5 * (dtaup[d][K*m+k] + dtaum[d][K*m+k]);
            }
        }
    }
}

void update_source_grid(const int D, const int I, const int K, const int M, 
                        double *E1, double *T6, double **source)
/**
 * \brief
 *
 * The function calculates the source function array.  Currently, the source 
 * function is just the specific intensity defined by Planck's Law because we
 * are assuming local thermodynamic equilibrium.
 *
 * \param[in] D       the total number of depth grid points
 * \param[in] I       the total number of energy and angular dependent grid
 *                    points
 * \param[in] K       the total number of energy  grid points
 * \param[in] M       the total number of angular dependent grid points
 * \param[in] E1      the energy grid (units: keV)
 * \param[in] T6      the temperature profile (units: 10^6 K)
 * \param[out] source the source function
 *                    (units: 10^22 erg cm^-2 s^-1 keV^-1 sr^-1)
 */                       
{
    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++)
        {
            for(int m = 0; m < M; m++) source[d][K*m+k] = planck_func(E1[k], 
                                                          T6[d]);
        }
    }
}

double planck_func(double E1, double T6)
/**
 * \brief
 *
 * This function calculates the value of the Planck spectrum (energy units). 
 * The return value has units of 10^22 erg cm^-2 keV^-1 s^-1 ster^-1.
 *
 * \param[in] E1 the energy value 
 * \param[in] T6 the temperature value
 * 
 * \return Outputs the value of Planck's formula at E1 and T6 as a double.
 */
{
    double out = 5.040385 * pow(E1, 3.0) / (exp(11.60445 * E1 / T6) - 1.0);

    return(out);
}

void calc_dBdtau_grid(const int D, const int K, const int M, 
                        double teff, double *dBdtau, double *E1, 
                        double *h, double *mu, double *T6, double **chi)
/**
 * \brief
 *
 * This function calulates the dBdtau term and all subsquent coefficients in
 * front and outputs an array of length I to be used for the iminus boundary 
 * condition at depth in rte().
 *
 * \param[in]  D      the total number of depth grid points
 * \param[in]  K      the total number of energy grid points
 * \param[in]  M      the total number of angular-dependent grid points
 * \param[in]  teff   the effective temperature of the NSA (units: 10^6 K)
 * \param[out] dBdtau the dBdtau term in BC of iminus at depth 
 *                    (units: 10^22 erg cm^-2 s^-1 keV^-1 sr^-1)
 * \param[in]  E1     the energy grid (units: keV)
 * \param[in]  h      the integration weights of the energy grid 
 *                    (units: keV)
 * \param[in]  mu     the angular dependent grid (units: None)
 * \param[in]  T6     the temperature profile (units: 10^6 K)
 * \param[in]  chi    the extinction coefficient (units: None)
 */
{
    double *dBdT = malloc_double_1d(K);
    double *f    = malloc_double_1d(K);

    for(int k = 0; k < K; k++)
    {
        double E14   = pow(E1[k], 4.0);
        double T6m2  = pow(T6[D-1], -2.0);
        double boltz = exp(11.60445 * E1[k] / T6[D-1]);
        double Sm2   = pow(boltz - 1, -2.0);

        dBdT[k] = 5.849089e-5 * E14 * T6m2 * boltz * Sm2;

        f[k]    = dBdT[k] / chi[D-1][k];
    }

    double ival = quadrature(K, f, h);
    double fbarc = 4.51245e-4 * pow(teff, 4.0);
    double dTdtauTh = 3.0 * fbarc / ival;
    
    for(int k = 0; k < K; k++)
    {
        for(int m = 0; m < M; m++) dBdtau[K*m+k] = 2.0 * mu[m] * dBdT[k] *
                                                   dTdtauTh / chi[D-1][k];
    }

    free_1d(dBdT);
    free_1d(f);
}
