/*
 * RaschModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef RASCHMODEL_H_
#define RASCHMODEL_H_

#include <model/parameter/ParameterModel.h>

class RaschModel: public ParameterModel {
public:
	// Constructor
	RaschModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void successProbability(DimensionModel *);

	// Getters and Setters
	map<Parameter, Matrix<double> *> getParameterSet() ;
	void setParameterSet(map<Parameter, Matrix<double> *>);
	double getProbability (int, int);

	// Destructor
	virtual ~RaschModel();
};

#endif /* RASCHMODEL_H_ */
