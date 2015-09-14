/* 
 * File:   extinction.h
 * Author: Justin Blythe
 *
 * Created on August 30, 2012, 9:53 PM
 */

void update_chi_grid(const int D, const int K, double *E1, double *rho, 
                     double *T6, double **chi);
double chi_ff(double E1, double rho, double T6);
double   g_ff(double gamma2,  double u);

