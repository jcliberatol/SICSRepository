/*
 * EM.cpp
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#include "EMEstimation.h"

EMEstimation::EMEstimation() {
	// TODO Auto-generated constructor stub

}

EMEstimation::~EMEstimation() {
	// TODO Auto-generated destructor stub
}
/*
 * Model is set.
 */
void EMEstimation::setModel(Model* Model){
	this->model=Model;
}
/*
 * Sets or selects the initial values
 */
//Sets
void EMEstimation::setInitialValues(map<Parameter, Matrix<double>*> parameterSet){
	model->getParameterModel()->setParameterSet(parameterSet);
}
//Selects
void EMEstimation::setInitialValues(string method){
	/*TODO
	 * Possible methods
	 * ANDRADE
	 * OSPINA
	 * RANDOM
	 *
	 * The default method is OSPINA
	 */
}

void EMEstimation::stepE () {
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

	PatternMatrix* data = dynamic_cast <PatternMatrix * > (model->getItemModel()->getDataset()) ;
	//Pattern iterator is data->iterator
	//Item number
	const double items = data->countItems();
	//Success probability matrix is obtained via pm->getProbability(int,int)
	ParameterModel* pm = model->getParameterModel();
	//Thetas
	Matrix<double>* thetas = model->getDimensionModel()->getLatentTraitSet()->getTheta();
	//Amount of nodes
	const int q = thetas->nC();
	//Weights
	Matrix<double>* weights = model->getDimensionModel()->getLatentTraitSet()->getWeight();
	//A Matrix
	Matrix<double>* A = model->getParameterModel()->getParameterSet()[a];
	//B Matrix
	Matrix<double>* B = model->getParameterModel()->getParameterSet()[b];
	//C Matrix
	Matrix<double>* C = model->getParameterModel()->getParameterSet()[c];
	//Auxiliar array for the nodes
	long double faux[q];
	long double sum = 0.0;
	int pat = 0; // N of patterns
	//Restart f and r to zero
	f->reset();
	r->reset();

	for (data->resetIterator();data->checkEnd();data->iterate()){

	}
	/*
	std::map<boost::dynamic_bitset<>, long int>::const_iterator iterator;

	long double faux[q];
	long double sum = 0.0;
	int pat = 0; // N of patterns

	memset(f, 0.0, q * sizeof(long double));
	memset(r, 0.0, q * items * sizeof(long double));

	for ( iterator = u.U.begin(); iterator != u.U.end(); ++iterator ) {

		pat++;
		//Initialize faux in 1 to later calculate the productory
		for ( int k = 0; k < q; k++ ) {
			faux[k] = 1;
		}

		//Calculate g*(k) for all the k's
		//first calculate the P for each k and store it in the array f aux
		for (int k = 0; k < q; k++) {

			//Calculate the p (iterate over the items in the productory)
			for ( unsigned int i = 0; i < items; i++ ) {
				faux[k] = faux[k]
						* successProbability(theta[k], a[i], b[i], c[i],
								(bool) iterator->first[i]);
			}

			//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
			//Now multiply by the weight

			faux[k] = faux[k] * nodeWeight[k];
		}

		//compute the total of the p*a' (denominator of the g*)
		sum = 0.0;
		for (int k = 0; k < q; k++) {
			sum += faux[k];
		}

		for (int k = 0; k < q; k++) {
			faux[k] = faux[k] / sum;	//This is g*_j_k

			//Multiply the f to the frequency of the pattern
			faux[k] = ((long double) iterator->second) * faux[k];
			f[k] = faux[k] + f[k];

			//Now selectively add the faux to the r
			for ( unsigned int i = 0; i < items; i++ ) {

				if (iterator->first[i]) {
					r_ki(k,i)= r_ki(k,i) + faux[k];

				} // if

			} // for

		} // for

	}*/

} //end E step

void EMEstimation::stepM(){
	/*
	 * TODO Step M for EML3M
	 */
}
void EMEstimation::estimate(){
	/*
	 * TODO Estimate
	 */
}


//Sets the method of getting boundary conditions and applying them
void EMEstimation::setBoundaryConditions(string option, string filename){

}
//Check if all the conditions are met for running the model, can report an error to a logger
void EMEstimation::checkRunningConditions(){

}
	//Sets the optimization algorithm
void EMEstimation::setOptimizationAlgorithm(string algorithm){

}
	//Sets the reporter for the trace
void EMEstimation::setTrace(string filename){

}
void EMEstimation::setTrace(Trace trace){

}


