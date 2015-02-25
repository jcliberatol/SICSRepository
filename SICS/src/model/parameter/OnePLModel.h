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
#include <type/Constant.h>
#include <cmath>

class OnePLModel : public ParameterModel {

public:
	// Constructor
	OnePLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	inline void successProbability(DimensionModel *, QuadratureNodes *);
	inline static double successProbability(double*, double*);
	static double logLikelihood(double*, double*, int, int);
	static void gradient (double* , double* , int , int , double* );
	// Getters and Setters
	double *** getParameterSet() ;
	void setParameterSet(double***);
	double getProbability (int, int);
	void getParameters(double * );
	void printParameterSet(ostream&);
	string getStringParameters();
	// Destructor
	virtual ~OnePLModel();
};


#endif /* ONEPLMODEL_H_ */
