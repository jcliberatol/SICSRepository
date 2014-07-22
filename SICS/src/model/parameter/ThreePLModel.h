/*
 * ThreePLModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef THREEPLMODEL_H_
#define THREEPLMODEL_H_

#include <typeinfo>
#include <model/parameter/ParameterModel.h>
#include <model/item/DichotomousModel.h>
#include <model/item/PolytomousModel.h>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultiUniDimModel.h>
#include <type/PatternMatrix.h>
#include <type/Constant.h>

class ThreePLModel: public ParameterModel {
	double successProbability(double, double, double, double);
public:
	// Constructor
	ThreePLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void successProbability(DimensionModel *);
	double LogLikelihood(double*, double*, int, int);

	// Getters and Setters
	map<Parameter, Matrix<double> *> getParameterSet() const;
	void setParameterSet(map<Parameter, Matrix<double> *>);
	double getProbability(int, int);

	// Destructor
	virtual ~ThreePLModel();
};

#endif /* THREEPLMODEL_H_ */
