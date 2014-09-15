/*
 * TwoPLModel.h
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#ifndef TWOPLMODEL_H_
#define TWOPLMODEL_H_

#include <model/parameter/ParameterModel.h>

class TwoPLModel: public ParameterModel {
public:
	// Constructor
	TwoPLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void successProbability(DimensionModel *);

	// Getters and Setters
	map<Parameter, Matrix<double> *> getParameterSet() ;
	void setParameterSet(map<Parameter, Matrix<double> *>);
	double getProbability (int, int);

	// Destructor
	virtual ~TwoPLModel();
};

#endif /* TWOPLMODEL_H_ */
