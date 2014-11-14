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
	estimator = NULL;
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


	//Discriminate by models
	if(Model->Modeltype()==Constant::THREE_PL){
	estimator = new EM3PL(); //Initializes estimator
	}

}
/**
 * Sets the initial values for the estimation, use this for inputting a matrix as initial values
 */
void EMEstimation::setInitialValues(map<Parameter, Matrix<double>*> parameterSet) {
	estimator->setInitialValues(parameterSet,model);
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
void EMEstimation::setInitialValues(int method) {
	estimator->setInitialValues(method,model);
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
	while (!convergenceSignal) {
		cout << "Iteration " << iterations << endl;
		estimator->stepE(model,f,r,quadNodes);
		estimator->stepM(model,f,r,quadNodes);
		convergenceSignal = model->itemParametersEstimated;
		iterations++;
		if (iterations > 200) {
			convergenceSignal = true;
			cout<<"more than 200 iters, stop"<<endl;
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
