/*
 * ThreePLModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/ThreePLModel.h>

ThreePLModel::ThreePLModel() {
	// TODO Auto-generated constructor stub

}

void ThreePLModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {
}

void ThreePLModel::setInitialPars(map<Parameter, Matrix<double> >* parameterSet) {
}

void ThreePLModel::calculateInitialPars() {
}

void ThreePLModel::successProbability() {
}

const map<Parameter, Matrix<double> >* ThreePLModel::getParameterSet() const {
	return parameterSet;
}

void ThreePLModel::setParameterSet(map<Parameter, Matrix<double> >* pair) {
	this->parameterSet = parameterSet;
}

ThreePLModel::~ThreePLModel() {
	// TODO Auto-generated destructor stub
}

