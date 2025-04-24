#ifndef MinFunc_hpp
#define MinFunc_hpp

typedef double (*MinimizeFxn) (double *);

double   LocalMin (double a, double b, double eps, double t, MinimizeFxn f, double *px);
double   PrAxis (double tol, double h, int nn, double xx[], MinimizeFxn f, double v[], double work[], int maxIterations);
void     BracketMinimum (double *pA, double *pB, MinimizeFxn func);
double   RandomNumber (long int *seed);

#endif
