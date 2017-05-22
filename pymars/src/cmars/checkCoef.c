#include "stdio.h"
#include <stdlib.h> 
#include <math.h>
#include "cmars.h"
#include "utilities.h"

checkCoefReturn checkCoef(double coef, double eps, int kr, double sw, double v, double we, double se, double dy, double *D, double *DY)
{
    /*This function checks to see if the proposed coefficient is best.*/
    checkCoefReturn cc = {0, 0.0};
    double dv, a, b;

    dv = v-pow(we,2)/sw;
    if (dv > 0.0)
    {
        a = inner(&D[1], &DY[1], kr+1);
        b = inner(&D[1], &D[1], kr+1);

        b = dv-b;
        if (b > eps*dv)
        {
            b = -pow(dy-a, 2)/b;
            /*this check is a very important one in the algorithm. if it
            passes the new hockey stick might be accepted and the response
            surface changes. In the original code there was no accounting for 
            numerical issues. The addition of 1.e-10 was added only after
            some painful debugging. The effect of adding this is to say 
            that the new coefficient (b) is taken as the new one if it's
            really smaller.*/
            if (b + 1.e-10< coef)
                {
                cc.newCoef = 1;
                cc.coef = b;
                return cc;  /*new coefficient*/
                }
        }
    }
    return cc;
}