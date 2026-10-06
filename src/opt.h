/*  opt.h
 *
 *  Optimizers used by ratePlacer: Brent's method in one dimension (Brent1D, via LocalMin) and
 *  in several dimensions (minimize_brent, via PrAxis), both from minfunc.c, plus a coarse
 *  golden-section search (GoldenSection_rough).
 */
#include "jph.h"
#include "minfunc.h"

#define INF DBL_MAX

#define TOLER_PASS_2        praxisTol
#define MAX_STEP_SIZE_2     1.0

double praxisTol = 0.00000001;

// Increase the tolerance value to make this a less stringent golden section to make faster for calls of golden section that don't need high accuracy of the likelihood
double GoldenSection_rough(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x,double z[]), double in_otherstuff[], int n_otherstuff)
{
	int i;
	double a = lowbound[1], b = upbound[1];
	double invphi = (sqrt(5) - 1) / 2;
	double invphi2 = (3 - sqrt(5)) / 2;
	double h = b - a;
	double c = a + invphi2 * h;
	double d = a + invphi * h;
	double yc, yd, tol = 0.000001; 	// same as in ameoba
	double *otherstuff;

	//first make one vector that contains bounds and other parameters
	otherstuff = malloc((2*n+n_otherstuff)*(sizeof(double)));
	for (i=0; i<n_otherstuff; i++)
	{
		otherstuff[i] = in_otherstuff[i];
		//printf("otherstuff %d %.16f\n", i, otherstuff[i]); 
	}
	for (i=0; i<n; i++)
	{
		otherstuff[2*i+n_otherstuff] = lowbound[i+1]; //the bounds start counting at 1
		otherstuff[2*i+n_otherstuff+1] = upbound[i+1]; //the bounds start counting at 1
							       //printf("Bounds: %lf %lf\n",otherstuff[2*i+n_otherstuff] ,otherstuff[2*i+n_otherstuff+1]);
	}

	yc = (*fun)(c, otherstuff);
	yd = (*fun)(d, otherstuff);

	while ( h > tol ){
		if (yc < yd) 
		{
			b = d;
			d = c;
			yd = yc;
			h = invphi * h;
			c = a + invphi2 * h;
			yc = (*fun)(c, otherstuff);
		}
		else 
		{
			a = c;
			c = d;
			yc = yd;
			h = invphi * h;
			d = a + invphi * h;
			yd = (*fun)(d, otherstuff);
		}
	}


	if (yc < yd)
	{
		newinvecter[1] = (a + d)/2;
		return(yc);
	}
	else
	{
		newinvecter[1] = (b + c)/2;
		return(yd);
	}
}

/*
 * Adapter for Brent1D: bridges GoldenSection-style callback
 * double (*fun)(double x, double z[]) to MinimizeFxn double (*)(double *, void *).
 */
typedef struct {
	double (*fun)(double x, double z[]);
	double *otherstuff;
} Brent1D_ctx;

static double Brent1D_bridge(double *x, void *extra_data)
{
	Brent1D_ctx *ctx = (Brent1D_ctx *)extra_data;
	return ctx->fun(*x, ctx->otherstuff);
}

/*
 * 1D Brent's method (LocalMin) with the same calling convention as GoldenSection.
 * Faster on smooth functions due to quadratic interpolation (superlinear convergence).
 * Drop-in replacement: result stored in newinvecter[1], returns minimum value.
 */
double Brent1D(double newinvecter[], double lowbound[], double upbound[], int n,
               double (*fun)(double x, double z[]), double in_otherstuff[], int n_otherstuff)
{
	int i;
	double px, result;

	/* Build the same combined otherstuff array layout as GoldenSection uses */
	double *otherstuff = malloc((2 * n + n_otherstuff) * sizeof(double));
	for (i = 0; i < n_otherstuff; i++)
		otherstuff[i] = in_otherstuff[i];
	for (i = 0; i < n; i++) {
		otherstuff[2 * i + n_otherstuff]     = lowbound[i + 1];
		otherstuff[2 * i + n_otherstuff + 1] = upbound[i + 1];
	}

	Brent1D_ctx ctx = {fun, otherstuff};
	result = LocalMin(lowbound[1], upbound[1], sqrt(DBL_EPSILON), 1e-8,
	                  Brent1D_bridge, &px, &ctx);

	newinvecter[1] = px;
	free(otherstuff);
	return result;
}

double minimize_brent(double newinvecter[], int n, double (*fun)(double x[], void*), int maxIterations, void *extra_data) {
    double minusLnL;

    double *directions = (double*)malloc(sizeof(double) * n * n);
    if (!directions)
    {
        printf ("Could not allocate directions (%lu)\n", sizeof(double) * n * n);
        exit (1);
    }
    
    double *powellWork = (double*)malloc(sizeof(double) * 6 * n);
    if (!powellWork)
    {
        printf ("Could not allocate powellWork (%lu)\n", sizeof(double) * 6 * n);
        exit (1);
    }

    minusLnL = PrAxis(TOLER_PASS_2, MAX_STEP_SIZE_2, n, newinvecter, *fun, directions, powellWork, maxIterations, extra_data);

    free(directions);
    free(powellWork);

    return(minusLnL);
}