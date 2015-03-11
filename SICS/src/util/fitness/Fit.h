/*
 * Fit.h
 *
 *  Created on: Mar 10, 2015
 *      Author: anmrdz
 */

#ifndef UTIL_FITNESS_FIT_H_
#define UTIL_FITNESS_FIT_H_
#include <type/LatentTraits.h>
#include <type/PatternMatrix.h>
#include <cmath>
#include <type/Constant.h>


void Fit(double**  LL,double** P,double** Q,LatentTraits* scores, Matrix<double> data,Model* model ){

	bool ** pattern_list = scores->pm->getBitsetList();
	int nitems = data.nC();
	int nscores = scores->pm->matrix.size();
	int ninds = data.nR();
	double scoresTot [ninds];
	int i, j, k;
	int npatt = 0;

	for(i = 0; i < ninds ; i++){
		for(j = 0; j < nscores; j++){
			npatt = 0;

			for(k = 0; k < nitems ; k++){
				if(data(i, k) == (pattern_list[j][k] > 0)? 1: 0)
					npatt++;
			}

			if(npatt == nitems){
				scoresTot[i] = (*(scores->traits))(j, 0);
				break;
			}
		}
	}


	double*** parameterSet = model->getParameterModel()->getParameterSet();
	int model_type = model->type;
	double a, b, c, cp;

	for(j = 0; j < ninds; j++){
		for(i = 0; i < nitems; i++){
			switch(model->type){
			case Constant::RASCH_A1 :
				a = 1;
				b = parameterSet[0][0][i];
				c = 0;
				break;
			case Constant::RASCH_A_CONSTANT :
				a = parameterSet[0][0][0];
				b = parameterSet[1][0][i];
				c = 0;
				break;
			case Constant::TWO_PL :
				a = parameterSet[0][0][i];
				b = parameterSet[1][0][i];
				c = 0;
				break;
			case Constant::THREE_PL:
				a = parameterSet[0][0][i];
				b = parameterSet[1][0][i];
				c = parameterSet[2][0][i];
				break;
			}

			cp = log(c/(1-c));
			P[j][i] = 1 / (exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-Constant::D_CONST*(a*scoresTot[j]+ b))));
			Q[j][i] = 1 - P[j][i];
		}
	}

	for(i = 0; i < nitems; i++){
		for(j = 0; j < ninds; j++){
			if(data(j, i) == 1){
				LL[j][i] = P[j][i];
			}else{
				LL[j][i] = Q[j][i];
			}
		}
	}
}

#endif /* UTIL_FITNESS_FIT_H_ */
