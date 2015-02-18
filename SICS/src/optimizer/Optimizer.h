/*
 * Optimizer.h
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_
#include <optimizer/BFGSOptimizer.h>

/**
 * Wrapper for the other optimizers, take general functions for the function to optimize, the vector gradient function and the matrix hessian function
 * The parameters are :
	 * double * args (Arguments over which the function optimizes)
	 * double * pars (Arguments in whose the function depends but are not optimized)
	 * int nargs Number of arguments
	 * int npars Number of parameters
	 * double * return (Return of the function is put in this array.)
 * */
class Optimizer {
public:
	void searchOptimal(double (*functionPtr)(double*,double*,int,int),
			void (*gradientPtr)(double*,double*,int,int,double*),
			void (*HessianPtr)(double*,double*,int,int,double*),
			double* args, double* pars, int nargs, int npars, double* hessiana);
	~Optimizer();

};

#endif /* OPTIMIZER_H_ */
