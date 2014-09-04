/*
 * EM.cpp
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#include "EMEstimation.h"

EMEstimation::EMEstimation() {
	iterations = 100;
	model = NULL;
	f = NULL;
	r = NULL;
	logger = NULL;
	optim = NULL;
	convergenceSignal = false;
	optim = new Optimizer();

}

EMEstimation::~EMEstimation() {
	if (f != NULL) {
		delete f;
	}
	if (r != NULL) {
		delete r;
	}
	if (logger != NULL) {
		delete logger;
	}
	if (optim != NULL) {
		delete optim;
	}
}
/*
 * Model is set.
 */
void EMEstimation::setModel(Model* Model) {
	int q;
	int It;
	this->model = Model;
	q = model->getDimensionModel()->getLatentTraitSet()->getTheta()->nC();
	It = model->getItemModel()->countItems();

	f = new Matrix<double>(1, q);
	r = new Matrix<double>(q, It);

}
/*
 * Sets or selects the initial values
 */
//Sets
void EMEstimation::setInitialValues(
		map<Parameter, Matrix<double>*> parameterSet) {
	model->getParameterModel()->setParameterSet(parameterSet);
}
//Selects
void EMEstimation::setInitialValues(string method) {
	/*TODO
	 * Possible methods
	 * ANDRADE
	 * OSPINA
	 * RANDOM
	 *
	 * The default method is OSPINA
	 */
}

void EMEstimation::stepE() {
	/*
	 * What we need
	 * q
	 * a pattern iterator
	 * item number
	 * success probability matrix
	 * thetas
	 * weights
	 * parameter set
	 */
	//Dataset by patterns
	PatternMatrix* data =
			dynamic_cast<PatternMatrix *>(model->getItemModel()->getDataset());
	//Pattern iterator is data->iterator
	//Item number
	const double items = data->countItems();
	//Success probability matrix is obtained via pm->getProbability(int,int)
	ParameterModel* pm = model->getParameterModel();
	//Thetas
	Matrix<double>* thetas =
			model->getDimensionModel()->getLatentTraitSet()->getTheta();
	//Amount of nodes
	const int q = thetas->nC();
	//Weights
	Matrix<double>* weights =
			model->getDimensionModel()->getLatentTraitSet()->getWeight();
	//A Matrix
	Matrix<double>* A = model->getParameterModel()->getParameterSet()[a];
	//B Matrix
	Matrix<double>* B = model->getParameterModel()->getParameterSet()[d];
	//C Matrix
	Matrix<double>* C = model->getParameterModel()->getParameterSet()[c];
	//Auxiliar array for the nodes
	long double faux[q];
	long double sum = 0.0;
	//Restart f and r to zero
	f->reset();
	r->reset();
	//cout<<(* data)<<endl;
	model->successProbability();
	//cout<<(*pm->probabilityMatrix)<<endl<<"pm"<<endl;

	for (data->resetIterator(); !data->checkEnd(); data->iterate()) {
		//Initialize faux in 1 to later calculate the productory
		for (int k = 0; k < q; k++) {
			faux[k] = 1;
		}
		//Calculate g*(k) for all the k's
		//first calculate the P for each k and store it in the array f aux
		for (int k = 0; k < q; k++) {
			//Calculate the p (iterate over the items in the productory)
			for (unsigned int i = 0; i < items; i++) {
				double prob = pm->getProbability(k, i);
				//cout<<data->getCurrentBitSet()[items-i-1];
				if (!data->getCurrentBitSet()[items - i - 1]) {
					prob = 1 - prob;
				}
				faux[k] = faux[k] * prob;
			}		//cout<<endl;
			//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
			//Now multiply by the weight
			faux[k] = faux[k] * (*weights)(0, k);
		}
		//cout<<(*weights);
		//compute the total of the p*a' (denominator of the g*)
		sum = 0.0;
		for (int k = 0; k < q; k++) {
			sum += faux[k];
			//cout<<faux[k]<<" ";
		}		//cout<<"Da sum ist : "<<sum<<endl;
		for (int k = 0; k < q; k++) {
			faux[k] = faux[k] / sum;	//This is g*_j_k
			//if(k==0)cout<<faux[k]<<"faux after div"<<endl;
			//Multiply the f to the frequency of the pattern
			faux[k] = ((long double) data->getCurrentFrequency()) * faux[k];
			//if(k==0)cout<<faux[k]<<"faux after freq"<<endl;
			(*f)(0, k) += faux[k];
			//if(k==0)cout<<(*f)(0,k)<<"The f "<<((long double) data->getCurrentFrequency())<<endl;
			//Now selectively add the faux to the r
			for (unsigned int i = 0; i < items; i++) {
				if (data->getCurrentBitSet()[items - i - 1]) {
					(*r)(k, i) = (*r)(k, i) + faux[k];
				} // if
			} // for
		} // for

	}
	//cout<<(*r)<<endl;
} //end E step

void EMEstimation::stepM() {
	/*
	 * TODO Step M for EML3M
	 */
	//Step M implementation using the BFGS Algorithm
	/*
	 * What we need
	 * fptr the pointer to loglik
	 * gprt the pointer to gradient
	 * hessptr the pointer to hessian matrix calculatrix
	 * args the a,b, and c
	 * pars, the other parameters q and stuff
	 * nargs, npars, sizes.
	 */
	//fptr
	double (*fptr)(double*, double*, int, int);
	void (*gptr)(double*, double*, int, int, double*);
	void (*hptr)(double*, double*, int, int, double*);
	fptr = &ThreePLModel::logLikelihood;
	gptr = &ThreePLModel::gradient;
	hptr = NULL;
	//cout<<"Address : "<<&gptr<<" "<<&hptr<<endl;
	int It = model->getItemModel()->getDataset()->countItems();
	int q = model->getDimensionModel()->getLatentTraitSet()->getTheta()->nC();
	double args[3 * It];
	double pars[2 + 2 * q + q * It];
	int nargs = 3 * It;
	int npars = 2 + 2 * q + q * It;
	//filling args
	int nA = 0;
	// Obtain a
	//A Matrix
	Matrix<double>* A = model->getParameterModel()->getParameterSet()[a];
	//B Matrix
	Matrix<double>* B = model->getParameterModel()->getParameterSet()[d];
	//C Matrix
	Matrix<double>* C = model->getParameterModel()->getParameterSet()[c];

	Matrix<double> DA(*A);
	Matrix<double> DB(*B);
	Matrix<double> DC(*C);

	for (int i = 0; i < It; i++) {
		args[nA] = (*A)(0, i);
		nA++;
	}

	// Obtain b
	for (int i = 0; i < It; i++) {
		args[nA] = (*B)(0, i);
		nA++;
	}

	// Obtain c
	for (int i = 0; i < It; i++) {
		args[nA] = (*C)(0, i);
		nA++;
	}
	//Filling pars
	int nP = 0;
	// Obtain q
	pars[nP] = q;
	nP++;
	// Obtain I
	pars[nP] = It;
	nP++;
	// Obtain theta
	//Thetas

	Matrix<double>* thetas =
			model->getDimensionModel()->getLatentTraitSet()->getTheta();
	for (int k = 0; k < q; k++) {
		pars[nP] = (*thetas)(0, k);	//TODO correct indexing on this and nearby matrices
		nP++;
	}
	// Obtain f
	for (int k = 0; k < q; k++) {
		pars[nP] = (*f)(0, k);
		nP++;
	}
	// Obtain r
	for (int k = 0; k < q; k++) {
		for (int i = 0; i < It; i++) {
			pars[nP] = (*r)(k, i);
			nP++;
		}
	}
	nargs = nA;
	npars = nP;
	/*
	 * Chooses the method
	 * method 1 is NR
	 * method 2 is BFGS
	 */
	int method = 1;

	if (method == 1) {
		//Newton Raphson
		double grad[3 * It];
		double hess[3 * 3 * It];
		ThreePLModel::gradient(args, pars, nargs, npars, grad);
		ThreePLModel::Hessian(args, pars, nargs, npars, hess);
		//TODO CHANGE TO LOGGER
		/*
		 cout<<"Gradient calculated"<<endl;
		 for (int p = 0; p < 3; ++p) {
		 for (int g = 0; g < It; ++g) {
		 cout<<grad[p*It+g]<<" ";
		 }cout<<endl;
		 }
		 cout<<"Hessian Calculated"<<endl;
		 cout<<endl;

		 for (int var = 0; var < It; ++var) {
		 for(int s = 0 ; s < 9 ; ++s){
		 cout<<hess[var*9+s]<<"  ";
		 }cout<<endl;
		 }cout<<endl;
		 */
		for (int i = 0; i < It; ++i) {
			double targs[3];
			//fill the args for the item
			targs[0] = args[i];
			targs[1] = args[It + i];
			targs[2] = args[2 * It + i];
			double tgrad[3];
			double thess[9];
			//Pass the item number through the tunnel in memory
			tgrad[0] = It;
			thess[0] = It;
			//Create the gradient pointer
			void (*tgptr)(double*, double*, int, int, double*);
			void (*thptr)(double*, double*, int, int, double*);
			tgptr = &ThreePLModel::itemgradient;
			thptr = &ThreePLModel::itemHessian;
			//optimize with these parameters changing the args
			newton(tgptr, thptr, targs, hess, grad, 3, i, 100, tgrad, thess);
			//cout<<"Made it this far";
			//update the args
			args[i] = targs[0];
			args[It + i] = targs[1];
			args[2 * It + i] = targs[2];

			//cout<<"Args : 	"<<args[i]<<" "<<args[It+i]<<" "<<args[2*It+i]<<" "<<endl;
		}
	}
	if (method == 2) {
		//BFGS
		optim->searchOptimal(fptr, gptr, hptr, args, pars, nargs, npars);
	}

	// Now pass the optimals to the Arrays.

	nA = 0;
	// Obtain a
	for (int i = 0; i < It; i++) {
		(*A)(0, i) = args[nA++];
		if (fabs((*A)(0, i)) > abs(10)) {
			(*C)(0, i) = 0.852;
			cout << "A reset." << endl;
		}

	}
	// Obtain b
	for (int i = 0; i < It; i++) {
		(*B)(0, i) = args[nA++];
		if (fabs((*B)(0, i)) > abs(-50)) {
			(*C)(0, i) = 0.5;
			cout << "B reset." << endl;
		}
	}
	// Obtain c
	for (int i = 0; i < It; i++) {
		(*C)(0, i) = args[nA++];
		if (fabs((*C)(0, i)) > abs(600)) {
			(*C)(0, i) = -1.7364;
			cout << "C reset." << endl;
		}
	}
	//Boundary regularize the arguments
	//C = -1.7346
	//B = 0.5;
	//A = 0.851

	//Obtain the deltas
	//Perform substracts
	double maxDelta = 0;
	double meanDelta = 0;
	int DeltaC = 0;
	for (int v1 = 0; v1 < It; ++v1) {
		DA(0, v1) = DA(0, v1) - (*A)(0, v1);
		DB(0, v1) = DB(0, v1) - (*B)(0, v1);
		DC(0, v1) = DC(0, v1) - (*C)(0, v1);
		meanDelta = +fabs(DA(0, v1));
		meanDelta = +fabs(DB(0, v1));
		meanDelta = +fabs(DC(0, v1));
		DeltaC += 3;
		if (fabs(DA(0, v1)) > maxDelta) {
			maxDelta = fabs(DA(0, v1));
		}
		if (fabs(DB(0, v1)) > maxDelta) {
			maxDelta = fabs(DB(0, v1));
		}
		if (fabs(DC(0, v1)) > maxDelta) {
			maxDelta = fabs(DC(0, v1));
		}
	}
	meanDelta = meanDelta / DeltaC;
	if (meanDelta < 0.0001) {
		convergenceSignal = true;
	}
	cout << "Max Delta : " << maxDelta << endl;
	cout << "Mean Delta : " << meanDelta << endl;
	cout << "MATS : " << endl << (*A) << (*B) << (*C) << endl;
	//And set the parameter sets
	map<Parameter, Matrix<double> *> parSet;
	parSet[a] = A;
	parSet[d] = B;
	parSet[c] = C;
	// llenar las tres matrices
	model->getParameterModel()->setParameterSet(parSet);

}

void EMEstimation::estimate() {
	/*
	 * TODO Estimate
	 */
	//Transform B and C
	cout << "Item impression" << endl;
	for (int i = 0; i < model->getItemModel()->countItems(); ++i) {
		double qa = (*model->getParameterModel()->getParameterSet()[a])(0, i);
		double qb = (*model->getParameterModel()->getParameterSet()[d])(0, i);
		double qc = (*model->getParameterModel()->getParameterSet()[c])(0, i);
		//(*model->getParameterModel()->getParameterSet()[d])(0,i)= -qb*qa;
		(*model->getParameterModel()->getParameterSet()[c])(0, i) = log(
				qc / (1 - qc));
	}
	iterations = 0;
	while (!convergenceSignal) {
		cout << "Iteration " << iterations << endl;
		stepE();
		//cout<<"____________________________________________________________________________"<<endl;
		//cout<<*model->getParameterModel()->getParameterSet()[a]
		//    <<*model->getParameterModel()->getParameterSet()[d]
		//    <<*model->getParameterModel()->getParameterSet()[c];
		stepM();
		iterations++;
		if (iterations > 100) {
			convergenceSignal = true;
		}
	}
	//Transform B
	//Transform C
	for (int i = 0; i < model->getItemModel()->countItems(); ++i) {
		double qa = (*model->getParameterModel()->getParameterSet()[a])(0, i);
		double qb = (*model->getParameterModel()->getParameterSet()[d])(0, i);
		double qc = (*model->getParameterModel()->getParameterSet()[c])(0, i);
		//(*model->getParameterModel()->getParameterSet()[d])(0,i)= -qb/qa;
		double ec = exp(qc);
		(*model->getParameterModel()->getParameterSet()[c])(0, i) = ec
				/ (1 + ec);
	}
	cout
			<< "______________FINAL VALUES FOR A  B AND THE C________________________"
			<< endl;
	cout << *model->getParameterModel()->getParameterSet()[a]
			<< *model->getParameterModel()->getParameterSet()[d]
			<< *model->getParameterModel()->getParameterSet()[c];
}
//Check if all the conditions are met for running the model, can report an error to a logger
void EMEstimation::checkRunningConditions() {

}
//Sets the optimization algorithm
void EMEstimation::setOptimizationAlgorithm(string algorithm) {

}
//Sets the reporter for the trace
void EMEstimation::setTrace(string filename) {

}
void EMEstimation::setTrace(Trace trace) {

}

int EMEstimation::getIterations() const {
	return (iterations);
}
