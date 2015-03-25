/*
 * ItemFit.h
 *
 *  Created on: 27/02/2015
 *      Author: anmrodriguezre
 */

#ifndef UTIL_FITNESS_ITEMFIT_H_
#define UTIL_FITNESS_ITEMFIT_H_
#include <type/LatentTraits.h>
#include <type/PatternMatrix.h>
#include <cmath>
#include <type/Constant.h>
#include <cstring>
#include <util/fitness/Fit.h>

void itemFit(LatentTraits* scores, Matrix<double> data,double*** parameterSet, int model_type, double * itemsf){

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

	Fit(LL, P, Q, scores, data, parameterSet, model_type);

	double sum = 0;

	for(j = 0; j < nitems; j++){
		sum = 0;
		for(i = 0; i < ninds ; i++){
			sum += log(LL[i][j]);
		}
		LL[i - 1][j] = sum;
	}

	double sigmaCuad[nitems], mu[nitems];
	memset(sigmaCuad, 0, sizeof(sigmaCuad));
	memset(mu, 0, sizeof(mu));

	for(i = 0; i < nitems; i++){
		for(j = 0; j < ninds; j++){
			mu[i] += log(P[j][i]) * P[j][i] + log(Q[j][i]) * Q[j][i];
			sigmaCuad[i] += P[j][i] * Q[j][i] * (log(P[j][i]/Q[j][i]) * log(P[j][i]/Q[j][i]));
		}
	}


	for(j = 0; j < nitems; j++){
		itemsf[j] = (LL[ninds - 1][j] - mu[j])/sqrt(sigmaCuad[j]);
		cout << itemsf[j] << "\n";
	}
}

#endif /* UTIL_FITNESS_ITEMFIT_H_ */
