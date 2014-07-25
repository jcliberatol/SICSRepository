/*
 * EM.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef EM_H_
#define EM_H_
#include <estimation/classical/ClassicalEstimation.h>
#include <string>
#include <trace/Trace.h>
#include <type/Matrix.h>
#include <type/PatternMatrix.h>
//#include <optimizer/BFGSOptimizer.h>

class EMEstimation : public ClassicalEstimation{
public:

	/*
	 * Order of execution
	 *
	 * Load the model
	 * Set or Select method of initial values
	 * Set the boundary conditions
	 * Set the optimization algorithm
	 * Sets a logger to report the results of the trace
	 * Checks the running conditions
	 * Estimates
	 */
	EMEstimation();
	virtual ~EMEstimation();
	//Executes the estimation
	void estimate();
	//Starts the E Step
	void stepE();
	//Starts the M Step
	void stepM();
	//Sets the model to estimate
	void setModel(Model* model);
	//Check if all the conditions are met for running the model, can report an error to a logger
	void checkRunningConditions();
	//Sets the optimization algorithm
	void setOptimizationAlgorithm(string algorithm);
	//Sets the reporter for the trace
	void setTrace(string filename);
	void setTrace(Trace trace);
	//Sets the initial values
	void setInitialValues(map<Parameter, Matrix<double>* > parameterSet);
	void setInitialValues(string method);

private:
	//F and R Matrices, remember to set to zero and open memory in process
	Matrix<double>* f;
	Matrix<double>* r;
	//Holds the trace for the logger outputs
	Trace* logger;
	//Model on which the algorithm operates
	Model* model;
	//Algorithm of optimization used
	Optimizer* optim;
	bool convergenceSignal;
};

#endif /* EM_H_ */


//
