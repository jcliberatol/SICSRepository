/*
 * EM.cpp
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#include <estimation/classical/EMEstimation.h>
#include <util/util.h>

EMEstimation::EMEstimation() {
	iterations = 0;
	model = NULL;
	f = NULL;
	r = NULL;
	logger = NULL;
	optim = NULL;
	convergenceSignal = false;
	optim = new Optimizer();
	quadNodes = NULL;

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
/**
 * Sets the model to be estimated, currently only supports 3PL model
 */
void EMEstimation::setModel(Model* Model) {
	int q;
	int It;
	this->model = Model;
	q = quadNodes->size();
	It = model->getItemModel()->countItems();

	f = new Matrix<double>(1, q);
	r = new Matrix<double>(q, It);

}
/**
 * Sets the initial values for the estimation, use this for inputting a matrix as initial values
 */
void EMEstimation::setInitialValues(map<Parameter, Matrix<double>*> parameterSet) {
	model->getParameterModel()->setParameterSet(parameterSet);
}
/**
 * Sets the initial values according to a method of calculating the values
 * Possible methods :
	 * ANDRADE,
	 * OSPINA,
	 * RANDOM,
	 *
	 * The default method is OSPINA , this is the fastest method according to the SICS calculations
 */
void EMEstimation::setInitialValues(string method) {
	/*TODO
	 * Possible methods
	 * ANDRADE
	 * OSPINA
	 * RANDOM
	 *
	 * The default method is OSPINA
	 */
	if(!method.compare("RANDOM")){
		std::srand(std::time(0)); // use current time as seed for random generator
		int items = model->getParameterModel()->getParameterSet()[a]->nC();
		for (int i = 0 ; i < items ; i++){
			(*model->getParameterModel()->getParameterSet()[a])(0, i)= randomd()*2;
			//fill b
			(*model->getParameterModel()->getParameterSet()[d])(0, i)= randomd()*4-2 ;
			//fill c
			int m = 4;
			(*model->getParameterModel()->getParameterSet()[c])(0, i)= randomd()*(2/(double)m);

			cout<<"i : "<<i<<" a "
								<<(*model->getParameterModel()->getParameterSet()[a])(0, i)<<" b "
								<<(*model->getParameterModel()->getParameterSet()[d])(0, i)<<" c "
								<<(*model->getParameterModel()->getParameterSet()[c])(0, i)<<endl;
		}
	}

	if(!method.compare("ANDRADE")){
	//Andrade method
	int items = model->getParameterModel()->getParameterSet()[a]->nC();
	//sums of the patterns
	int totalscores = 0 ;
	int *itemscores = new int [items];
	memset(itemscores,0,sizeof(int)*items);
	double *covariances = new double [items];
	memset(covariances,0,sizeof(double)*items);
	double variance = 0;
	PatternMatrix* data = dynamic_cast<PatternMatrix *>(model->getItemModel()->getDataset());
	double Ni = (double)data->countIndividuals();
	for (data->resetIterator(); !data->checkEnd(); data->iterate()) {
		double df = (double)data->getCurrentFrequency();
		double bs = (double)data->getCurrentBitSet().count();
		for (int i = 0 ; i < items ; i++){
			if(data->getCurrentBitSet()[i]){
				itemscores[i]+=df;
			}
		}
		totalscores += bs*df;
	}
	cout<<"Score total : "<<totalscores<<endl;
	//calculate variances and covariances
	for (data->resetIterator(); !data->checkEnd(); data->iterate()){
		double df = (double)data->getCurrentFrequency();
		double bs = (double)data->getCurrentBitSet().count();
		for (int i = 0 ; i < items ; i++){
			if(data->getCurrentBitSet()[i]){
				covariances[i]+=((1-itemscores[i]/Ni)*(1-bs/items))*df;
			}
		}
		variance+=((bs-((double)totalscores/Ni))*(bs-((double)totalscores/Ni)))*df;
	}
	variance /= Ni;
	for (int i = 0 ; i < items ; i++){
		covariances[i] /= (Ni-1);
		cout<<"cov : "<<i<<" "<<covariances[i]<<endl;
	}
	//Now calculate the standard deviations for the sums and the items
	long double*stddevs = new long double [items];
	memset(stddevs,0,sizeof(long double)*items);
	long double*pearson = new long double [items];
	memset(pearson,0,sizeof(long double)*items);
	long double*pis = new long double [items];
	memset(pis,0,sizeof(long double)*items);
	for (int i = 0 ; i < items ; i++){
			double avg= totalscores/Ni;
			stddevs[i]=stdDev_bin(itemscores[i],Ni,avg);
			pis[i]=itemscores[i]/Ni;
			pearson[i]=(covariances[i]/(stddevs[i]*std::sqrt(variance)));
			//fill a sqrt(pCoef * pCoef / (1.0 - pCoef * pCoef));
			(*model->getParameterModel()->getParameterSet()[a])(0, i)= std::sqrt((pearson[i]*pearson[i])/(1/pearson[i]*pearson[i]));
			//fill b
			(*model->getParameterModel()->getParameterSet()[d])(0, i)=normalInverse(pis[i]);
			//fill c
			int m = 4;
			(*model->getParameterModel()->getParameterSet()[c])(0, i)= 1 / (double)m;//TODO CHANGE BY CONSTANT FROM CONST.H FILE
			cout<<"i : "<<i<<" a "
					<<(*model->getParameterModel()->getParameterSet()[a])(0, i)<<" b "
					<<(*model->getParameterModel()->getParameterSet()[d])(0, i)<<" c "
					<<(*model->getParameterModel()->getParameterSet()[c])(0, i)<<endl;
	}
	}
}

/**
 * Step E of the EM method, this step takes the actual
 * estimation of the parameters, and produces a f and a r
 * matrices used in the maximization step
 *
 * TODO : PARALLELIZABLE FOR
 */
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
	PatternMatrix* data = dynamic_cast<PatternMatrix *>(model->getItemModel()->getDataset());
	//Pattern iterator is data->iterator
	//Item number
	const double items = data->countItems();
	//Success probability matrix is obtained via pm->getProbability(int,int)
	ParameterModel* pm = model->getParameterModel();
	//Thetas
	Matrix<double>* thetas = quadNodes->getTheta();
	//Amount of nodes
	const int q = quadNodes->size();
	//Weights
	Matrix<double>* weights =quadNodes->getWeight();
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
	//Calculates the success probability for all the nodes.
	model->successProbability(quadNodes);

	//TODO CAREFULLY PARALLELIZE FOR
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
				if (!data->getCurrentBitSet()[items - i - 1]) {
					prob = 1 - prob;
				}
				faux[k] = faux[k] * prob;
			}
			//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
			//Now multiply by the weight
			faux[k] = faux[k] * (*weights)(0, k);
		}
		//compute the total of the p*a' (denominator of the g*)
		sum = 0.0;
		for (int k = 0; k < q; k++) {
			sum += faux[k];
		}
		for (int k = 0; k < q; k++) {
			faux[k] = faux[k] / sum;	//This is g*_j_k
			//Multiply the f to the frequency of the pattern
			faux[k] = ((long double) data->getCurrentFrequency()) * faux[k];

			(*f)(0, k) += faux[k];
			//Now selectively add the faux to the r
			for (unsigned int i = 0; i < items; i++) {
				if (data->getCurrentBitSet()[items - i - 1]) {
					(*r)(k, i) = (*r)(k, i) + faux[k];
				} // if
			} // for
		} // for

	}
} //end E step

/**
 * Executes the maximization step using the inputted or default optimizer
 * currently only supporting BFGS and newton algorithms
 * */
void EMEstimation::stepM() {
	/*
	 */
	//Step M implementation using the BFGS Algorithm
	//Function pointers to represent the loglikelihood, gradient and hessian
	double (*fptr)(double*, double*, int, int);
	void (*gptr)(double*, double*, int, int, double*);
	void (*hptr)(double*, double*, int, int, double*);
	fptr = &ThreePLModel::logLikelihood;
	gptr = &ThreePLModel::gradient;
	hptr = NULL;
	int It = model->getItemModel()->getDataset()->countItems();
	int q = quadNodes->size();
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
	//Thetas

	Matrix<double>* thetas =quadNodes->getTheta();
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
	//BFGS
		optim->searchOptimal(fptr, gptr, hptr, args, pars, nargs, npars);

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
	if (meanDelta < 0.0001 and maxDelta < 0.001) {
		convergenceSignal = true;
	}
	//And set the parameter sets
	map<Parameter, Matrix<double> *> parSet;
	parSet[a] = A;
	parSet[d] = B;
	parSet[c] = C;
	// llenar las tres matrices
	model->getParameterModel()->setParameterSet(parSet);

}

/**
 * Main loop of the EM estimation
 * orchestrates the parameters for the estimation, and holds the estimation
 * for the iterations.
 *
 * TODO : read maxiterations as a input parameter , idea : calculate the max iterations depending on the items
 * TODO : Output last estimation onto the json for recovery in the program.
 */
void EMEstimation::estimate() {
	/*
	 * TODO Estimate
	 */
	//Transform B and C
	for (int i = 0; i < model->getItemModel()->countItems(); ++i) {
		double qa = (*model->getParameterModel()->getParameterSet()[a])(0, i);
		double qb = (*model->getParameterModel()->getParameterSet()[d])(0, i);
		double qc = (*model->getParameterModel()->getParameterSet()[c])(0, i);
		(*model->getParameterModel()->getParameterSet()[c])(0, i) = log(
				qc / (1 - qc));
	}
	iterations = 0;
	EMEstimator* estimator = new EM3PL();
	while (!convergenceSignal) {
		cout << "Iteration " << iterations << endl;

		//stepE();
		estimator->stepE(model,f,r,quadNodes);
		//stepM();
		estimator->stepM(model,f,r,quadNodes);
		convergenceSignal = model->itemParametersEstimated;
		iterations++;
		if (iterations > 200) {
			convergenceSignal = true;
			cout<<"more than 200 iters, stopp"<<endl;
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
			<< "___________________CONVERGENCE VALUES________________________"
			<< endl;
	   cout << *model->getParameterModel()->getParameterSet()[a]
			<< *model->getParameterModel()->getParameterSet()[d]
			<< *model->getParameterModel()->getParameterSet()[c]
			<<"______________________________________________________________"<<endl;
}
/**Check if all the conditions are met for running the model, can report an error to a logger*/
void EMEstimation::checkRunningConditions() {

}
/**Sets the optimization algorithm*/
void EMEstimation::setOptimizationAlgorithm(string algorithm) {

}
/**Sets the reporter for the trace*/
void EMEstimation::setTrace(string filename) {

}
void EMEstimation::setTrace(Trace trace) {

}
/**Returns the iterations that took the estimation to obtain an answer*/
int EMEstimation::getIterations() const {
	return (iterations);
}

QuadratureNodes* EMEstimation::getQuadratureNodes() const {
	return (quadNodes);
}

void EMEstimation::setQuadratureNodes(QuadratureNodes* nodes) {
	this->quadNodes = nodes;
}
