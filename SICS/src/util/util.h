/*
 * util.h
 *
 *  Created on: 26/06/2014
 *      Author: jliberato
 */

#ifndef UTIL_H_
#define UTIL_H_
#include <util/asa111.hpp>

/**
 * Functions part one
 * Includes the next functions
 * StdDev_bin  : Calculates the standard deviation for a binary vector.
 * Determinant3_3 : Calculates the determinant of a 3 by 3 , 2d matrix.
 * logTransform : Transforms from the 0 to 1 scale to the real scale (read logit)
 * ExpTransform : Transforms from the real scale to the 0 , 1 scale (read antilogit)
 * PCheck : checks if a number is a probability, if it is one returns a very big probability, if its zero, a very small but not null.
 */

inline double randomd(){

	int random_variable = std::rand();
	return ((double)((double)random_variable/(double)RAND_MAX));
}

inline long double stdDev_bin( int tsum, int tN){
	long double avg;
	long double N = (long double)tN;
	long double sum = (long double)tsum;
	avg = sum/N;
	return (std::sqrt((sum*(1 - avg)*(1 - avg)+(N-sum)*(avg*avg))/N));
}

inline long double stdDev_bin( int tsum, int tN, double avg){
	long double N = (long double)tN;
	long double sum = (long double)tsum;
	return (sqrt((((1-avg)*(1-avg)*sum)+((-avg)*(-avg)*(N-sum)))/N));
}

inline double normalInverse(double point){
	int err = 0;
	return(ppnd(point,&err));
}

//   v1 := actual parameters
//   v1 := parameters of 1 step back
//   v2 := parameters of 2 step back
inline void ramsay(double *v0, double *v1, double *v2, int size)
{
	double dX[size], dX2[size], d2X2[size], accel, accelI, numerator = 0.0, denominator = 0.0;
	for ( int i = 0; i < size; i++ )
	{
		dX[i] = v0[i]  - v1[i];
		dX2[i] = v1[i] - v2[i];
		d2X2[i] = dX[i] - dX2[i];
	}
	for ( int i = 0; i < size; i++ )
	{
		numerator += dX[i]*dX[i];
		denominator += d2X2[i]*d2X2[i];
	}
	accel = 1 - sqrt(numerator/denominator);
	if ( accel < -5.0) accel = -5;
	accelI = 1-accel;
	for ( int i = 0; i < size; i++ ) v0[i] = accelI*v0[i]+accel*v1[i];
}


#endif /* UTIL_H_ */
