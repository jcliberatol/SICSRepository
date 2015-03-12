/*
 * util.h
 *
 *  Created on: 26/06/2014
 *      Author: jliberato
 */

#ifndef UTIL_H_
#define UTIL_H_
#include <util/asa111.hpp>

#define last_it (*args_hist)[2]
#define prev_it (*args_hist)[1]
#define b_prev_it (*args_hist)[0]

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
	double dX[size];
	double dX2[size];
	double d2X2[size];
	double accel;
	double numerator = 0.0;
	double denominator = 0.0;

	for (int i = 0; i < size; i++) {
		dX[i] = last_it[i] - prev_it[i];
		dX2[i] = prev_it[i] - b_prev_it[i];
		d2X2[i] = dX[i] - dX2[i];

		numerator += dX[i] * dX[i];
		denominator += d2X2[i] * d2X2[i];
	}

	accel = 1 - sqrt(numerator / denominator);

	if (accel < -5.0)
		accel = -5;

	for (int i = 0; i < size; i++)
		last_it[i] = (1 - accel) * last_it[i] + accel * prev_it[i];
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
