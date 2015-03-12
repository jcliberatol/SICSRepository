/*
 * PersonFit.h
 *
 *  Created on: Mar 11, 2015
 *      Author: anmrodriguezre
 */

#ifndef UTIL_FITNESS_PERSONFIT_H_
#define UTIL_FITNESS_PERSONFIT_H_

#include <type/LatentTraits.h>
#include <type/PatternMatrix.h>
#include <cmath>
#include <type/Constant.h>
#include <cstring>
#include <util/fitness/Fit.h>

void personFit(LatentTraits* scores, Matrix<double> data,Model* model ){

	int nitems = data.nC();
	int ninds = data.nR();
	int i, j;
	double* LL[ninds];
	double* P[ninds];
	double* Q[ninds];

	for(i = 0; i < ninds; i++){
		LL[i] = new double[nitems];
		P[i] = new double[nitems];
		Q[i] = new double[nitems];
	}

	Fit(LL, P, Q, scores, data, model);

	double sum = 0;

	for(i = 0; i < ninds ; i++){
		sum = 0;
		for(j = 0; j < nitems; j++){
			sum += log(LL[i][j]);
		}
		LL[i][j - 1] = sum;
	}

	double sigmaCuad[ninds], mu[ninds];
	memset(sigmaCuad, 0, sizeof(sigmaCuad));
	memset(mu, 0, sizeof(mu));

	for(i = 0; i < nitems; i++){
		for(j = 0; j < ninds; j++){
			mu[j] += log(P[j][i]) * P[j][i] + log(Q[j][i]) * Q[j][i];
			sigmaCuad[j] += P[j][i] * Q[j][i] * (log(P[j][i]/Q[j][i]) * log(P[j][i]/Q[j][i]));
		}
	}

	double Z3[nitems];

	for(i = 0; i < ninds; i++){
		Z3[i] = (LL[i][nitems - 1] - mu[i])/sqrt(sigmaCuad[i]);
		cout << Z3[i] << "\n";
	}
}

#endif /* UTIL_FITNESS_PERSONFIT_H_ */
