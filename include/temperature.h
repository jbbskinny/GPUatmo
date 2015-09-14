/* 
 * File:   temperature.h
 * Author: Justin Blythe
 *
 * Created on September 5, 2012, 10:43 AM
 */

void unsold_lucy_temp_correct(const int D, const int I, const int K,
                              const int M, const int maxiteration,
                              const double g, const double teff, 
                              double *E1, double *h, double *mu, double *rho,
                              double *T6, double *tauTh, double *w, 
                              double **tau, double **taum, double **taup,
                              double **chi, double **iminus, double **iplus, 
                              double **source);
void rosseland_temp_profile(const int D, const int K, const int maxiteration,
                            const double g, const double teff, double *E1, 
                            double *h, double *rho, double *T6, double *tauTh);
void update_rossld_profile(const int D, const double teff,
                           double *T6, double *tau);
void calc_rossld_extinc(const int D, const int K, double *chiR, double *E1,
                        double *h, double *rho, double *T6);
void calc_rossld_depth(const int D, double *chiR, double *tauR, double *tauTh);

