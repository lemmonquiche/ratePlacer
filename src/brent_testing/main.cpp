#include <cmath>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "minfunc.h"
#include "jph.h"

#    define TOLER_PASS_1        0.001
#    define TOLER_PASS_2        praxisTol
#    define MAX_STEP_SIZE_1    1.0
#    define MAX_STEP_SIZE_2    1.0

double maximizeLike(void);
double CalcLGivenParamVect(double* x);

double *paramVector, *directions, brentEps, brentT, praxisTol = 0.00000001,
       *powellWork, *bestParamVector, bestLnL, minusLnL;
std::vector<double> data;

int main(void) {

	//RandomVariable& rng = RandomVariable::randomVariableInstance();

	//double mu = 3.0;
	//double sigma = 0.2;
	//int n = 1000;

	//double sum = 0.0;
	//for (int i=0; i<n; i++)
	//{
	//	//double x = Probability::Normal::rv(&rng, mu, sigma);
	//	double y = Probability::Uniform::rv(&rng, -100.0, 100.0);
	//	double z = Probability::Uniform::rv(&rng, -100.0, 100.0);
	//	double x = pow(z + 2*y -7, 2) + pow(2 * z + y - 5, 2);
	//	sum += x;
	//	data.push_back(x);
	//}
	//double mean = sum / data.size();
	//sum = 0.0;
	//for (int i=0; i<n; i++)
	//	sum += ((data[i] - mean) * (data[i] - mean));
	//double var = sum / (data.size()-1);

	//std::cout << "mean = " << mean << std::endl;
	//std::cout << "var = " << sqrt(var) << std::endl;

	double vals[] = { 0.1, 1.0 };
	std::cout << "test = " << CalcLGivenParamVect(vals) << std::endl;
	double lnL = maximizeLike();
	printf("lnL = %lf\n", lnL);

	//for (int i = 0; i<n; i++)
	//	std::cout << data[i] << ",";
	//std::cout << std::endl;

	return 0;
}

double maximizeLike(void) {

	int np = 2;

	paramVector = (double*)malloc(sizeof(double) * np);
	if (!paramVector)
	{
		printf ("Could not allocate paramVector (%lu)\n", sizeof(double) * np);
		exit (1);
	}
	//bestParamVector = (double*)malloc(sizeof(double) * np);
	//if (!bestParamVector)
	//{
	//	printf ("Could not allocate bestParamVector (%lu)\n", sizeof(double) * np);
	//	exit (1);
	//}
	directions = (double*)malloc(sizeof(double) * np * np);
	if (!directions)
	{
		printf ("Could not allocate directions (%lu)\n", sizeof(double) * np * np);
		exit (1);
	}
	powellWork = (double*)malloc(sizeof(double) * 6 * np);
	if (!powellWork)
	{
		printf ("Could not allocate powellWork (%lu)\n", sizeof(double) * 6 * np);
		exit (1);
	}

	//bestLnL = -1000000000;

	paramVector[0] = -500.0;
	paramVector[1] = 500.0;

	/* initial approximation */
	/*minusLnL = PrAxis(TOLER_PASS_1, MAX_STEP_SIZE_1, np, paramVector, CalcLGivenParamVect,
	  directions, powellWork);*/

	/* final estimation */
	minusLnL = PrAxis(TOLER_PASS_2, MAX_STEP_SIZE_2, np, paramVector, CalcLGivenParamVect, directions, powellWork);

	std::cout << "Best params: x: " << paramVector[0] << ", y: " << paramVector[1] << " with likelihood " << minusLnL << std::endl;

	std::cout << "Best params: x: " << bestParamVector[0] << ", y: " << bestParamVector[1] << " with likelihood " << minusLnL << std::endl;

	free (paramVector);
	free (bestParamVector);
	free (directions);
	free (powellWork);

	return minusLnL;
}

double CalcLGivenParamVect(double* x) {
	//if(x[0] < -100.0)
	//	x[0] = -100;
	//i
	std::cout << "CalcLGivenParamVect " << x[0] << " " << x[1] << " ";

	double lnL = pow(x[0] + 2*x[1] -7, 2) + pow(2 * x[0] + x[1] - 5, 2);

	std::cout << lnL << std::endl;

	return lnL;
}

//double CalcLGivenParamVect(double* x) {
//    double mu = x[0];
//    if (x[1] < 0.00000001)
//        x[1] = 0.00000001;
//    double sigma = x[1];
//    std::cout << "CalcLGivenParamVect " << mu << " " << sigma << " ";
//
//    double lnL = 0.0;
//    for (int i=0; i<data.size(); i++) 
//        lnL += Probability::Normal::lnPdf(mu, sigma, data[i]);
//
//    std::cout << lnL << std::endl;
//    
//    return -lnL;
//}

