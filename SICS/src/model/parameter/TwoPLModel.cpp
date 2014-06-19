/*
 * TwoPLModel.cpp
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/TwoPLModel.h>

TwoPLModel::TwoPLModel() {
	// TODO Auto-generated constructor stub

}

void TwoPLModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {
}

void TwoPLModel::setInitialPars(map<Parameter, Matrix<double> >* pair) {
}

void TwoPLModel::calculateInitialPars() {
}

void TwoPLModel::successProbability() {
}

const map<Parameter, Matrix<double> >* TwoPLModel::getParameterSet() const {
	return parameterSet;
}

void TwoPLModel::setParameterSet(map<Parameter, Matrix<double> >* parameterSet) {
	this->parameterSet = parameterSet;
}

TwoPLModel::~TwoPLModel() {
	// TODO Auto-generated destructor stub
}

