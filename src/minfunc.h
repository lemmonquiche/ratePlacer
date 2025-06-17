#ifndef MinFunc_hpp
#define MinFunc_hpp

typedef double (*MinimizeFxn) (double *, void *);

double   LocalMin (double a, double b, double eps, double t, MinimizeFxn f, double *px, void *extra_data);
double   PrAxis (double tol, double h, int nn, double xx[], MinimizeFxn f, double v[], double work[], int maxIterations, void *extra_data);
void     BracketMinimum (double *pA, double *pB, MinimizeFxn func, void *extra_data);
double   RandomNumber (long int *seed);

#endif
