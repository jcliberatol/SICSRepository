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

void EMEstimation::stepE(){
	/*
	 * TODO Step E for EML3M
	 */

}

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


