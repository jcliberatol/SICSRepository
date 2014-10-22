/*
 * TwoPLModel.h
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#ifndef TWOPLMODEL_H_
#define TWOPLMODEL_H_

#include <model/parameter/ParameterModel.h>
#include <type/Constant.h>
#include <cmath>
#include <typeinfo>
#include <model/dimension/UnidimensionalModel.h>

class TwoPLModel: public ParameterModel {
public:
	// Constructor
	TwoPLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void successProbability(DimensionModel *, QuadratureNodes *);
	static double successProbability(double, double, double, double);
	static double logLikelihood(double*, double*, int, int);

	// Getters and Setters
	map<Parameter, Matrix<double> *> getParameterSet() ;
	void setParameterSet(map<Parameter, Matrix<double> *>);
	double getProbability (int, int);

	// Destructor
	virtual ~TwoPLModel();
};

#endif /* TWOPLMODEL_H_ */
