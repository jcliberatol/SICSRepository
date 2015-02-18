#include <optimizer/Optimizer.h>

	void Optimizer::searchOptimal(double (*functionPtr)(double*,double*,int,int),
			void (*gradientPtr)(double*,double*,int,int,double*),
			void (*HessianPtr)(double*,double*,int,int,double*),
			double* args, double* pars, int nargs, int npars, double* hessiana){
	//TODO REAL OPTIMIZER
	//FOR NOW ONLY OPTIMIZE USING BFGS
	int r = bfgs(functionPtr,gradientPtr,args,pars,nargs,npars,1000, hessiana);
}

Optimizer::~Optimizer() {
	// TODO Auto-generated destructor stub
}
