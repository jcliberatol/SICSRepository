/*
 * EMEstimator.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 */

#ifndef EMESTIMATOR_H_
#define EMESTIMATOR_H_
#include <string>
#include <trace/Trace.h>
#include <type/Matrix.h>
#include <type/PatternMatrix.h>
#include <type/QuadratureNodes.h>
#include <optimizer/Optimizer.h>
#include <util/util.h>
//#include <estimation/classical/EMEstimators/EM3PL.h>

class EMEstimator {
public:
	EMEstimator(){}
	//Transforms the parameters before starting an estimation process
	virtual void transform(Model* m) = 0;
	//Transforms back the parameters after estimating them
	virtual void untransform(Model* m) = 0;
	//Step E needs the model , the f and r, and the thetas, besides from the data.
	virtual void stepE(Model* m, Matrix<double>* f, Matrix<double>* r,  QuadratureNodes* nodes) = 0;
	//Step M also needs the model, quad nodes, f and r
	virtual void stepM(Model* m, Matrix<double>* f, Matrix<double>* r,  QuadratureNodes* nodes) = 0 ;
	//in this cases the model are needed to be filled
	virtual void setInitialValues(int Method, Model*) = 0;
	virtual void setInitialValues(map<Parameter, Matrix<double>*> parameterSet, Model* m) = 0;
	virtual ~EMEstimator(){}
};

#endif /* EMESTIMATOR_H_ */
