/*
 * File: utilities.h
 * Author: Justin Blythe
 * 
 * Created November 28, 2012, 10:15 PM
 */

double residual(const int D, const double relerr, const double abserr,
                double *old, double *new);
void copy_array(const int size, double *copy, double *array);
double ** malloc_double_2d(const int dim1, const int dim2);
double * malloc_double_1d(const int dim1);
void free_2d(const int dim1, double **array);
void free_1d(double  *array);

