#include <optimizer/Optimizer.h>

	void Optimizer::searchOptimal(double (*functionPtr)(double*,double*,int,int),
			void (*gradientPtr)(double*,double*,int,int,double*),
			void (*HessianPtr)(double*,double*,int,int,double*),
			double* args, double* pars, int nargs, int npars){
	//TODO REAL OPTIMIZER
	//FOR NOW ONLY OPTIMIZE USING BFGS
	bfgs(functionPtr,gradientPtr,args,pars,nargs,npars,100);
}

Optimizer::~Optimizer() {
	// TODO Auto-generated destructor stub
}
