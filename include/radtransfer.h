/* 
 * File:   radtransfer.h
 * Author: Justin Blythe
 *
 * Created on October 5, 2012, 12:43 AM
 */

void solve_rte(const int D, const int I, const int K, const int M, double teff,
               double *E1, double *h, double *mu, double *T6, double **chi,
               double **dtau, double **dtaum, double **dtaup, double **iminus,
               double **iplus, double **source);
void update_dtau_grids(const int D, const int K, const int M, double *mu, 
                       double *tauTh, double **chi, double **dtau, 
                       double **dtaum, double **dtaup);
void update_source_grid(const int D, const int I, const int K, const int M, 
                        double *E1, double *T6, double **source);
double planck_func(double E1, double T6);
void calc_dBdtau_grid(const int D, const int K, const int M,
                      double teff, double *dBdtau, double *E1,
                      double *h, double *mu, double *T6, double **chi);
