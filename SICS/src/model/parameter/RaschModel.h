/*
 * RaschModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef RASCHMODEL_H_
#define RASCHMODEL_H_

#include <model/parameter/ParameterModel.h>
#include <typeinfo>
#include <model/item/DichotomousModel.h>
#include <model/item/PolytomousModel.h>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultiUniDimModel.h>
#include <type/PatternMatrix.h>
#include <type/Constant.h>
#include <cmath>
#include <type/QuadratureNodes.h>

class RaschModel: public ParameterModel {
private:
	QuadratureNodes* nodes;
public:

	static double successProbability(double, double);

	// Constructor
	RaschModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void successProbability(DimensionModel *, QuadratureNodes *);
	static double logLikelihood(double*, double*, int, int);
	void setEstimationNodes(QuadratureNodes*);
	static void gradient (double* , double* , int , int , double* );
	static void Ngradient(double* , double* , int , int , double* );
	static void NHessian(double* , double* , int , int , double* );

	// Getters and Setters
	double *** getParameterSet() ;
	void setParameterSet(double ***);
	double getProbability (int, int);
	string getStringParameters();

	// Destructor
	virtual ~RaschModel();
};

#endif /* RASCHMODEL_H_ */
