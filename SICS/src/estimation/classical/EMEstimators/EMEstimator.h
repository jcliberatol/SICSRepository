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
#include <type/Constant.h>

class EMEstimator {
public:
	Trace* profiler = 0;
	EMEstimator(){}
	EMEstimator(Model* m, QuadratureNodes* nodes, Matrix<double>* f, Matrix<double>* r){}
	//Transforms the parameters before starting an estimation process
	virtual void transform() = 0;
	//Transforms back the parameters after estimating them
	virtual void untransform() = 0;
	//Step E needs the model , the f and r, and the thetas, besides from the data.
	virtual void stepE() = 0;
	//Step M also needs the model, quad nodes, f and r
	virtual void stepM() = 0 ;
	//in this cases the model are needed to be filled
	virtual void setInitialValues(int , Model*) = 0;
	virtual void setInitialValues(double*** , Model*) = 0;
	void setProfiler(Trace* t){
		profiler = t;
	}
	virtual ~EMEstimator(){}
};

#endif /* EMESTIMATOR_H_ */
