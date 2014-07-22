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

void TwoPLModel::successProbability(DimensionModel *dimensionModel) {
}

map<Parameter, Matrix<double> *> TwoPLModel::getParameterSet() const {
	return (this->parameterSet);
}

void TwoPLModel::setParameterSet(map<Parameter, Matrix<double> *> pair) {
	this->parameterSet = parameterSet;
}

double TwoPLModel::getProbability(int node, int item) {
}

TwoPLModel::~TwoPLModel() {
	// TODO Auto-generated destructor stub
}

