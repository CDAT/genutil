/* This include file contains all of the structures & prototypes
 * associated with the mars functions.
 *
 */
typedef struct
{
    double v;
    double u;
    double t;
    double we;
    double se;
    double sy;
    double dy;
    double *ST;
    double *SU;
} computeCoefReturn;
computeCoefReturn computeCoef(double, double *, double *, double, double *, double,
                              long *, int, double, double *, int,  double *, int, int);

typedef struct
{
    int newCoef;
    double coef;
} checkCoefReturn;
checkCoefReturn checkCoef(double, double, int, double, double, double, double, double, double *, double *);
