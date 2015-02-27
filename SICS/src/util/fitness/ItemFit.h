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

void itemFit(LatentTraits* scores, PatternMatrix* data,Model* model ){
	int nitems = data->countItems();
	int nscores = (scores)->traits->nR();
	int ninds = data->countIndividuals();
	double scoresTot [nscores];
	int i, j, k;
	int npatt = 0;

	for(i = 0; i < ninds ; i++){
		for(j = 0; j < nscores; j++){
			npatt = 0;
			for(k = 0; k < nitems ; k++){
				if(data->bitset_list[i][k] == scores->pm->bitset_list[j][k])
					npatt++;
			}
			if(npatt == nitems){
				scoresTot[i] = (*(scores->traits))(j, 0);
			}
		}
	}

	double  LL[ninds][nitems], P[ninds][nitems], Q[ninds][nitems];
	double*** parameterSet = model->parameterModel->parameterSet;
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
			//TODO check qlogis function
			cp = log(c/(1-c));
			P[j][i] = 1 / (exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-Constant::D_CONST*(a*scoresTot[j]+ b))));
			Q[j][i] = 1 - P[j][i];
		}
	}

//	LL[data == 1] = P[data == 1]
//	  LL[data == 0] = Q[data == 0]
//	  LL = colSums(log(LL))
//
//
//	  mu = sigmaCuad = rep(0,nitems)
//	  for( i in 1:nitems){
//	    Pi = cbind(P[,i],Q[,i])
//	    logPi = log(Pi)
//	    mu[i] = sum(Pi * logPi)
//	    #sigmaCuad = sigmaCuad + Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2)
//	    sigmaCuad[i] = sum(Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2))
//
//	  }
//	  print("mu")
//	  print(mu)
//	  print(sigmaCuad)
//	  Z3 = (LL - mu) / sqrt(sigmaCuad)

}


#endif /* UTIL_FITNESS_ITEMFIT_H_ */
