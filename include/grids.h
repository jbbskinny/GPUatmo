/* 
 * File:   grids.h
 * Author: Justin Blythe
 *
 * Created on August 30, 2012, 9:55 PM
 */

void setup_tauTh_grid(const int decades, const int ppd, double init, 
                      double *tauTh);
void setup_E1_grid(const int decades, const int ppd, double init, double *E1,
                   double *h);
void setup_mu_grid(const int M, double *mu, double *w);
void update_density_profile(const int D, const double g, double *rho, 
                            double *T6, double *tauTh);
