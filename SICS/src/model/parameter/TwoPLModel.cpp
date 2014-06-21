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

const map<Parameter, Matrix<double> *>& TwoPLModel::getParameterSet() const {
	return this->parameterSet;
}

void TwoPLModel::setParameterSet(const map<Parameter, Matrix<double> *>& pair) {
	this->parameterSet = parameterSet;
}

TwoPLModel::~TwoPLModel() {
	// TODO Auto-generated destructor stub
}

