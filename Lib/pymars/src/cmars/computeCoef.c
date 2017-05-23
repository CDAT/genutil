/*
 * computeCoef.c
 *
 *  Created on: Dec 15, 2011
 *      Author: mcenerney1
 */
 #include <stdio.h>
 #include <stdlib.h> 
 #include <math.h>
 
#include "cmars.h"
#include "utilities.h"

computeCoefReturn computeCoef(double trialKnot, double *X, double *y, double yb, double *w, double sw, 
                              long *MM, int MM_len, double xt, double *partialEval, int kr, double *db, int db_nrows, int db_ncols)
{

    /*This function computes the local calculations involved in deciding
    a new knot and coefficient. It is called from bestKnot and once complete
    the accumulation of the local contributions is done before it is decided
    if the best one is found.*/
    

    /*printf("computeCoef in C\n");*/
    double h, xx, xd, su, st, yc, sj;
    int mj_ind, mj, i, j;
    double *pntr;
        
    computeCoefReturn cc = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    cc.ST = (double *)calloc(kr+1, sizeof(double));
    cc.SU = (double *)calloc(kr+1, sizeof(double));

    arrayInit(cc.ST, kr+1, 0.0);
    arrayInit(cc.SU, kr+1, 0.0);

    /*ST = numpy.zeros(shape = kr+1, dtype = FLOAT_DTYPE)
    SU = numpy.zeros(shape = kr+1, dtype = FLOAT_DTYPE)*/

    for (mj_ind=0; mj_ind<MM_len; mj_ind++)
    {
        mj = MM[mj_ind];
        h = partialEval[mj];
        if (w[mj] > 0.0 && h > 0.0)
        {
            xx = X[mj];
            xd = xx - trialKnot;
            su = w[mj]*h;
            st = su*xd;
            yc = y[mj] - yb;
            cc.dy = cc.dy + st*yc;
            cc.sy = cc.sy + su*yc;
            cc.we = cc.we + st;
            cc.se = cc.se + su;
            sj = w[mj]*h*h;
            cc.v = cc.v + sj*xd*xd;
            cc.t = cc.t + sj;
            cc.u = cc.u + sj*(xx - xt);
            /*ST[1:kr+1] = ST[1:kr+1] + st*db[mj,1:kr+1]
            SU[1:kr+1] = SU[1:kr+1] + su*db[mj,1:kr+1]
            
            printf("computeCoef in C, kr=%u\n", kr);
            fflush(NULL);*/
            pntr = &db[mj];
            /*linearUpdate(st, pntr+1, cc.ST+1, kr);
            linearUpdate(su, pntr+1, cc.SU+1, kr);
            linearUpdate(st, db, db_nrows, db_ncols, mj, cc.ST, 1, kr);
            linearUpdate(su, db, db_nrows, db_ncols, mj, cc.SU, 1, kr);*/
            
            printf("computeCoef in C, db_ncols=%u\n", db_ncols);
            fflush(NULL);   
            for (i=1; i<=kr; i++)
            {
                j = mj*db_ncols + i;
                cc.ST[i] += st*db[j];
                cc.SU[i] += su*db[j];
            }

        }
    }

    return cc;     /*return v, t, u, we, se, sy, dy, ST, SU*/
}