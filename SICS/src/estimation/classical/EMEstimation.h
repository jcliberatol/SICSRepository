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
#include <optimizer/Optimizer.h>
#include <trace/Trace.h>
#include <type/Matrix.h>
#include <type/PatternMatrix.h>

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
	//Sets the method of getting boundary conditions and applying them
	void setBoundaryConditions(string option, string filename);
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
	//Boundary matrix
	Matrix<double> *boundaryConditions;
	Optimizer* optim;
	//Algorithm of optimization used
};

#endif /* EM_H_ */


//
