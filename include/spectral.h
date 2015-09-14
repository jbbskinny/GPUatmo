/**
 * File: spectral.h
 * Author: Justin Blythe
 *
 * Created on March 11, 2013, 8:42 PM
 */
 
void update_jnu_grid(const int D, const int K, const int M, double *w,
                     double **iminus, double **iplus, double **jnu);
void update_j_profile(const int D, const int K, double *h, double *j,
                      double **jnu);
void update_chiJ_profile(const int D, const int K, double *chiJ, double *h,
                         double *j, double **chi, double **jnu);
void update_fnu_grid(const int D, const int K, const int M, double *mu, 
                     double *w, double **fnu, double **iminus, double **iplus);
void update_fbar_profile(const int D, const int K, double *fbar, double *h, 
                         double **fnu);
void update_chiF_profile(const int D, const int K, double *chiF, double *f,
                         double *h, double **chi, double **fnu);
void update_b_profile(const int D, double *b, double *T6);
void update_chiP_profile(const int D, const int K, double *b, double *chiP,
                         double *h, double **chi, double **source);
