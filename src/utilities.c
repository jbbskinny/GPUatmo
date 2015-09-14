/*
 * File: utilities.c
 * Author: Justin Blythe
 *
 * Created November 28, 2012, 10:05 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*****************************************************************************
 ***************************** Utilities *************************************
 *****************************************************************************/

double residual(const int N, const double relerr, const double abserr,
                double *old, double *new)
/**
 * \brief
 *
 * This function calculates the residual of an old and a new array.  Given the
 * total number of points, the relative error, and the absolute error, it 
 * calculates the maximum residual from the correct value.
 *
 * \param[in] N      the number of elements in array
 * \param[in] relerr the relative between old and new arrays
 * \param[in] abserr the absolute error
 * \param[in] old    the old array before update procedure was applied, length
 *                   N
 * \param[in] new    the new array, length N
 *
 * \return the maximum residual
 */                
{
    double dummy = 0;
    double resid;

    for(int n = 0; n < N; n++)
    {
        resid = (new[n] - old[n]) / (abserr + (relerr * new[n]));

        if(fabs(resid) > fabs(dummy)) dummy = resid;
    }
    
    return(dummy);
}

void copy_array(const int N, double *copy, double *array)
/**
 * \brief
 *
 * This function simply copies an array.
 *
 * \param[in]  N     the number of elements in the array
 * \param[out] copy  the array to be copied to
 * \param[in]  array the array to be copied
 */
{
    for(int n = 0; n < N; n++) copy[n] = array[n];
}

double ** malloc_double_2d(const int dim1, const int dim2)
/**
 * \brief
 * 
 * This function allocates a two dimensional array of doubles and returns a
 * pointer to the array.
 *
 * \param[in] dim1 the number of elements in the first dimension
 * \param[in] dim2 the number of elements in the second dimension
 *
 * \return the pointer to the two dimensional array
 *
 * \warning If memory is exceeded this function will quit the program, 
 *          returning the file name and line number.
 */
{
    double **array = (double**)malloc((size_t)sizeof(double) * dim1);

    if(array == NULL) 
    {
        printf("ERROR in %s at line %d.\n", __FILE__, __LINE__);
        exit(-1);
    }

    for(int i = 0; i < dim1; i++)
    {
        array[i] = (double*)malloc((size_t)sizeof(double) * dim2);

        if(array[i] == NULL)
        {
            printf("ERROR in %s at line %d.\n", __FILE__, __LINE__);
            exit(-1);
        }
    }

    return(array);
}

double * malloc_double_1d(const int dim1)
/**
 * \brief
 * 
 * This function allocates a one dimensional array of doubles and returns a
 * pointer to the array.
 *
 * \param[in] dim1 the number of elements in the first dimension
 *
 * \return the pointer to the one dimensional array of doubles
 *
 * \warning If memory is exceeded this function will quit the program, 
 *          returning the file name and line number.
 */
{
    double *array = (double*)malloc((size_t)sizeof(double) * dim1);

    if(array == NULL) 
    {
        printf("ERROR in %s at line %d.\n", __FILE__, __LINE__);
        exit(-1);
    }

    return (array);
}

void free_2d(const int dim1, double **array)
/**
 * \brief
 *
 * This function frees a two dimensional array.
 *
 * \param[in] dim1 the number of 1d array within the 2d array
 * \param[in] array the pointer to the 2d array
 */
{
    for(int i = 0; i < dim1; i++) free(array[i]);

    free(array);
}

void free_1d(double *array)
/**
 * \brief
 *
 * This function frees a one dimensional array.
 *
 * \param[in] array the pointer to the 1d array
 */
{
    free(array);
}

