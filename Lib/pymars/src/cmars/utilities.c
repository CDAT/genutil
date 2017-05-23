/*
 * utilities.c
 *
 *  Created on: Dec 21, 2011
 *      Author: mcenerney1
 */
#include "cmars.h"

void arrayInit(double *array, int length, double value)
{
    int i;
    for(i=0; i<length; i++)
        array[i] = value;
    return;
}

void linearUpdate(double slope, double *x, int x_nrows, int x_ncols, int row, double *y, int col0, int ncols)
{
    int i, j;
    for (i=col0; i<ncols; i++)
    {
        j = row*x_ncols + i;
        y[i] += slope*x[j];
    }
    return;
}

double inner(double x[], double y[], long length)
{
    long i;
    double sum = 0;
    for (i=0; i<length; i++)
        sum += x[i]*y[i];
    return sum;
}