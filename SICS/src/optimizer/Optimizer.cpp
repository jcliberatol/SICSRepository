#include <optimizer/Optimizer.h>

	void Optimizer::searchOptimal(double (*functionPtr)(double*,double*,int,int),
			void (*gradientPtr)(double*,double*,int,int,double*),
			void (*HessianPtr)(double*,double*,int,int,double*),
			double* args, double* pars, int nargs, int npars){
	//TODO REAL OPTIMIZER
	//FOR NOW ONLY OPTIMIZE USING BFGS
	int r = bfgs(functionPtr,gradientPtr,args,pars,nargs,npars,15);
	r  = 4;
	switch(r){
		case 1: cout<<"MAX ITER REACHED"<<endl; break;
		case 2: cout<<"BAD INIT VALUES"<<endl; break;
		case 3: cout<<"N_CONVERGENCE"<<endl; break;
		case 0:	cout<<"OPTIM SUCCESS"<<endl; break;
		case 4: break;
		default : cout<<"check the optimizer, an unexpected error is ocurring "<<endl; break;
	}
}

Optimizer::~Optimizer() {
	// TODO Auto-generated destructor stub
}
