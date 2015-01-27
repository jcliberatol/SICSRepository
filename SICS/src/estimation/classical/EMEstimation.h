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
#include <type/QuadratureNodes.h>
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <estimation/classical/EMEstimators/EM3PL.h>
#include <estimation/classical/EMEstimators/EM2PL.h>
#include <estimation/classical/EMEstimators/EM1PL.h>
#include <estimation/classical/EMEstimators/EM1PLAC.h>
//#include <optimizer/BFGSOptimizer.h>

/**
 *Classical estimation through the EM algorithm, generic for the models however must be called with specific
 *model object, the optimization algorithm can be any from the optimizers class.
 *EM is a iterative Estimation Maximization algorithm, please check the literature on how it works.
 *EM estimation requires a quadrature for the implementation of the integrals in the expectation step
 *this quadratures can be obtained from R, or using the supplied ones from the SICS binary quadratures
 *ranging from 1 quadrature node to 101 quadrature nodes.
 */
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
	void setInitialValues(double*** parameterSet);
	void setInitialValues(int method);
	int getIterations() const;
	QuadratureNodes* getQuadratureNodes() const;
	void setQuadratureNodes(QuadratureNodes* latentTraitSet);
	void setProfiler(Trace* t);
private:
	Trace* profiler;
	//F and R Matrices, remember to set to zero and open memory in process
	Matrix<double>* f;
	Matrix<double>* r;
	//Holds the trace for the logger outputs
	Trace* logger;
	//Estimator object to call estimation functions
	EMEstimator* estimator;
	//Model on which the algorithm operates
	Model* model;
	//Algorithm of optimization used
	Optimizer* optim;
	// true when convergence is reached
	bool convergenceSignal;
	//Number of iterations
	int iterations;
	//Sets used for the estimation
	QuadratureNodes *quadNodes;
};

#endif /* EM_H_ */
