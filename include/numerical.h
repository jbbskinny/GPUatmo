/* 
 * File:   numerical.c
 * Author: Justin Blythe
 *
 * Created on September 10, 2012, 1:19 PM
 */

double clenshaw(const int N, double x, double *C);
void setup_gauss_legendre(const int L, const double init, const double final, 
                          double *absys, double *weight);
double quadrature(const int L, double *f, double *weight);

