/*
 * Optimizer.h
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

class Optimizer {
public:
	virtual void Optimize();
	virtual ~Optimizer();
	//The optimizers general functions
	//The parameters are :
	/*
	 * double * args (Arguments over which the function optimizes)
	 * double * pars (Arguments in whose the function depends but are not optimized)
	 * int nargs Number of arguments
	 * int npars Number of parameters
	 * double * return (Return of the function is put in this array.)
	 * int Size of the return array
	 */
	void (*functionPtr)(double*,double*,int,int,double*,int);
	void (*gradientPtr)(double*,double*,int,int,double*,int);
	void (*HessianPtr)(double*,double*,int,int,double*,int);
};

#endif /* OPTIMIZER_H_ */
