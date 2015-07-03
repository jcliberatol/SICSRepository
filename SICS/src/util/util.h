/*
 * util.h
 *
 *  Created on: 26/06/2014
 *      Author: jliberato
 */

#ifndef UTIL_H_
#define UTIL_H_
#include <util/asa111.h>

#define iteration_0 (*args_hist)[2]
#define iteration_1 (*args_hist)[1]
#define iteration_2 (*args_hist)[0]

/**
 * Functions part one
 * Includes the next functions
 * StdDev_bin  : Calculates the standard deviation for a binary vector.
 * Determinant3_3 : Calculates the determinant of a 3 by 3 , 2d matrix.
 * logTransform : Transforms from the 0 to 1 scale to the real scale (read logit)
 * ExpTransform : Transforms from the real scale to the 0 , 1 scale (read antilogit)
 * PCheck : checks if a number is a probability, if it is one returns a very big probability, if its zero, a very small but not null.
 */

inline double randomd() {

	int random_variable = std::rand();
	return ((double) ((double) random_variable / (double) RAND_MAX));
}

inline long double stdDev_bin(int tsum, int tN) {
	long double avg;
	long double N = (long double) tN;
	long double sum = (long double) tsum;
	avg = sum / N;
	return (std::sqrt(
			(sum * (1 - avg) * (1 - avg) + (N - sum) * (avg * avg)) / N));
}

inline long double stdDev_bin(int tsum, int tN, double avg) {
	long double N = (long double) tN;
	long double sum = (long double) tsum;
	return (sqrt(
			(((1 - avg) * (1 - avg) * sum) + ((-avg) * (-avg) * (N - sum))) / N));
}

inline double normalInverse(double point) {
	int err = 0;
	return (ppnd(point, &err));
}

inline void ramsay(double *** args_hist, int size) {
#if DEBUG
	cout<<"Start Init values"<<endl;
#endif
	double dX[size];
	double dX2[size];
	double d2X2[size];
	double accel;
	double numerator = 0.0;
	double denominator = 0.0;
#if DEBUG
	cout<<"End Init values"<<endl;
#endif

#if DEBUG
	cout<<"Start for cicle"<<endl;
#endif

	for (int i = 0; i < size; i++) {
		dX[i] = iteration_0[i] - iteration_1[i];
		dX2[i] = iteration_1[i] - iteration_2[i];
		d2X2[i] = dX[i] - dX2[i];

		numerator += dX[i] * dX[i];
		denominator += d2X2[i] * d2X2[i];
	}
#if DEBUG
	cout<<"End for cicle"<<endl;
#endif
	accel = 1 - sqrt(numerator / denominator);

	if (accel < -5.0)
		accel = -5;

#if DEBUG
	cout<<"Start update values"<<endl;
#endif

	for (int i = 0; i < size; i++)
		iteration_0[i] = (1 - accel) * iteration_0[i] + accel * iteration_1[i];
#if DEBUG
	cout<<"End update values"<<endl;
#endif
}

inline void transformHessiana( double * inputHessiana, double ** outputHessiana, int size)
{
    for ( int i = 0; i < size; i++ )
    {
    	for ( int j = 0; j <= i; j++ )
    	{
    		outputHessiana[i][j] = inputHessiana[i+j];
    		outputHessiana[j][i] = inputHessiana[i+j];
    	}
    }
}

#endif /* UTIL_H_ */
