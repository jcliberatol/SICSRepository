/*
 * OnePLModel.h
 *
 *  Created on: Nov 16, 2014
 *      Author: anmrodriguezre
 */

#ifndef ONEPLMODEL_H_
#define ONEPLMODEL_H_

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

class OnePLModel : public ParameterModel {

private:
	QuadratureNodes* nodes;

public:

	static double successProbability(double, double);

	// Constructor
	OnePLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void successProbability(DimensionModel *, QuadratureNodes *);
	static double logLikelihood(double*, double*, int, int);
	void setEstimationNodes(QuadratureNodes*);
	static void gradient (double* , double* , int , int , double* );
	// Getters and Setters
	double *** getParameterSet() ;
	void setParameterSet(double***);
	double getProbability (int, int);
	void printParameterSet(ostream&);
	// Destructor
	virtual ~OnePLModel();
};


#endif /* ONEPLMODEL_H_ */
