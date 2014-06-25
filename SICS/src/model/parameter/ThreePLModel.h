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

class ThreePLModel: public ParameterModel {
public:
	// Constructor
	ThreePLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void setInitialPars(map<Parameter, Matrix<double> > *);
	void calculateInitialPars();
	void successProbability();

	// Getters and Setters
	const map<Parameter, Matrix<double> *>& getParameterSet() const;
	void setParameterSet(const map<Parameter, Matrix<double> *>&);

	// Destructor
	virtual ~ThreePLModel();
};

#endif /* THREEPLMODEL_H_ */
