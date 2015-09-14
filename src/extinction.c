/*
 * File: extinction.c
 * Author: Justin Blythe
 *
 * Created August 30, 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "extinction.h"
#include "numerical.h"

/*****************************************************************************
 ************************* Extinction: Free-free *****************************
 *****************************************************************************/

void update_chi_grid(const int D, const int K, double *E1, double *rho, 
                     double *T6, double **chi)
/**
 * \brief
 *
 * This function will initialize / update the chi grid. Where chi is a function
 * of depth and frequency. Note that the real extinction coefficient is in
 * units of Thompson extinction, thus it is unitless.
 *
 * \param[in]  D   the total number of depth dependent grid points
 * \param[in]  K   the total number of energy dependent grid points
 * \param[in]  E1  the energy grid (units: keV)
 * \param[in]  rho the density profile (units: g cm^-3)
 * \param[in]  T6  the temperature profile (units: 10^6 K)
 * \param[out] chi the extinction coefficient grid (units: None)
 */
{
    for(int d = 0; d < D; d++)
    {
        for(int k = 0; k < K; k++)
        {
            chi[d][k] = chi_ff(E1[k], rho[d], T6[d]);
        }
    }
}

double chi_ff(double E1, double rho, double T6)
/**
 * \brief
 *
 * This function calculates the free-free absorption coefficient, in units of 
 * Thompson extinction, thus it is unitless.  It is dependent on the energy, 
 * density, and temperature.  The only contribution to the extinction along the
 * beam direction is due to free-free absorption because the atmosphere of the
 * neutron star is completely ionized.
 *
 * \param[in] E1  the energy of a photon in the beam (units: keV)
 * \param[in] rho the density of the medium at a particular depth grid point
 *                (units: g cm^-3)
 * \param[in] T6  the temperature of the medium at a particular depth grid 
 *                point (units: 10^6 K)
 *
 * \return the extinction coefficient (units: none)
 */
{
    double u         = 11.60445 * (E1 / T6);
    double gamma2    = 0.15789 / T6;
    double prefactor = 1.86788;
    
    double chi = prefactor * g_ff(gamma2, u) * (1.0 - exp(-u)) * pow(rho, 1.0) 
                 / (pow(E1, 3.0) * pow(T6, 0.5));

    return(chi);
}

double g_ff(double gamma2, double u)  
/**
 * \brief
 *
 * This function calculates the gaunt factor, used to calculate the extinction 
 * coefficient. It uses the coefficient maxtrix of Chebyshev polynomials
 * determined by Hummer and then applies Clenshaw recurrance to evaluate the
 * function.
 *
 * \param[in] gamma2 an input that is the ionization divided by kT
 * \param[in] u      an input that is the photon energy divided by kT
 *
 * \return the calulated gaunt factor
 */
{
    const int N=8, M=11;
    
    double D[11][8] = 
    { 
        { 8.986940175e+0, -4.009515855e+0,  8.808871266e-1,  2.640245111e-2, 
           -4.580645915e-2, -3.568055702e-3,  2.827798067e-3,  3.365860195e-4},
        {-8.006936989e-1,  9.466021705e-1,  9.043402532e-2, -9.608451450e-2, 
           -1.885629865e-2,  1.050313890e-2,  2.800889961e-3, -1.078209202e-3},
        {-3.781305103e-1,  1.102726332e-1, -1.543619180e-2,  8.310561114e-3, 
            2.179620525e-2,  4.259726289e-3, -4.181588794e-3, -1.770208330e-3},
        { 1.877213132e-2, -1.004885705e-1, -5.483366378e-2, -4.520154409e-3, 
            8.366530426e-3,  3.700273930e-3,  6.889320423e-4,  9.460313195e-5},
        { 7.300158392e-2,  3.576785497e-3, -4.545307025e-3, -1.017965604e-2,
           -9.530211924e-3, -3.450186162e-3,  1.040482914e-3,  1.407073544e-3},
        {-1.744671550e-3,  2.864013856e-2,  1.903394837e-2,  7.091074494e-3,
           -9.668371391e-4, -2.999107465e-3, -1.820642230e-3, -3.874082085e-4},
        {-1.707268366e-2, -4.694254776e-3,  1.311691517e-3,  5.316703136e-3,
            5.178193095e-3,  2.451228935e-3, -2.277321615e-5, -8.182359057e-4},
        { 2.567331664e-4, -9.155339970e-3, -6.997479192e-3, -3.571518641e-3, 
           -2.096101038e-4,  1.553822487e-3,  1.509584686e-3,  6.212627837e-4},
        { 4.098322531e-3,  1.635218463e-3, -5.918883504e-4, -2.333091048e-3,
           -2.484138313e-3, -1.359996060e-3, -5.371426147e-5,  5.553549563e-4},
        { 3.837562402e-5,  2.938325230e-3,  2.393747064e-3,  1.328839809e-3, 
            9.135013312e-5, -7.137252303e-4, -7.656848158e-4, -3.504683798e-4},
        {-8.491991820e-4, -3.615327726e-4,  3.148015257e-4,  8.909207650e-4, 
            9.869737522e-4,  6.134671184e-4,  1.068883394e-4, -2.046080100e-4}
    };
    
    double C[8];
    double DD[11];
    
    double au = -1.9091; double bu = 0.09091;
    double ag = -1;      double bg = 1;

    double xg = log10(gamma2) / 3;
    double xu = (2 * log10(u) - 2.5) / 5.5;

    double yg = ((2.0 * xg) - ag - bg) / (bg - ag);
    double yu = ((2.0 * xu) - au - bu) / (bu - au);
    
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < M; i++) DD[i] = D[i][j];

        C[j] = clenshaw(M, yg, DD);
    }
    
    double out = clenshaw(N, yu, C);

    return (out);
}

