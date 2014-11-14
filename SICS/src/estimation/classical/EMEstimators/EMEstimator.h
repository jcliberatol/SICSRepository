/*
 * EMEstimator.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 */

#ifndef EMESTIMATOR_H_
#define EMESTIMATOR_H_
#include <estimation/classical/ClassicalEstimation.h>
#include <string>
#include <trace/Trace.h>
#include <type/Matrix.h>
#include <type/PatternMatrix.h>
#include <type/QuadratureNodes.h>
#include <optimizer/Optimizer.h>

class EMEstimator {
public:
	EMEstimator();
	//Step E needs the model , the f and r, and the thetas, besides from the data.
	void stepE(Model* m, Matrix<double>* f, Matrix<double>* r,  QuadratureNodes* nodes);
	//Step M also needs the model, quad nodes, f and r
	void stepM(Model* m, Matrix<double>* f, Matrix<double>* r,  QuadratureNodes* nodes);
	//in this cases the model are needed to be filled
	void setInitialValues(string Method);
	void setInitialValues(map<Parameter, Matrix<double>*> parameterSet);
	virtual ~EMEstimator();
};

#endif /* EMESTIMATOR_H_ */
