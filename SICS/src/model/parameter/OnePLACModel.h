#ifndef ONEPLACMODEL_H_
#define ONEPLACMODEL_H_

#include <typeinfo>
#include <model/parameter/ParameterModel.h>
#include <model/item/DichotomousModel.h>
#include <model/item/PolytomousModel.h>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultiUniDimModel.h>
#include <type/PatternMatrix.h>
#include <type/Constant.h>
#include <cmath>
#include <type/QuadratureNodes.h>

class OnePLACModel: public ParameterModel
{

private:
	QuadratureNodes* nodes;
public:

	static double successProbability(double, double, double);
	double successProbability(double, double *);
	// Constructor
	OnePLACModel();

	// Methods
	void setEstimationNodes(QuadratureNodes *);
	void successProbability(DimensionModel *, QuadratureNodes *);
	static double logLikelihood(double*, double*, int, int);
	static void gradient(double*, double*, int, int, double*);
	virtual void transform() {}
	virtual void untransform() {}
	// Getters and Setters
	double *** getParameterSet();
	void getParameters(double *);
	void setParameters(double *);
	void setParameterSet(double ***);
	double getProbability(int, int);
	void printParameterSet(ostream&);
};

#endif
