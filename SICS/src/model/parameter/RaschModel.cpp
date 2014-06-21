/*
 * RaschModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/RaschModel.h>

RaschModel::RaschModel() {
	// TODO Auto-generated constructor stub

}

void RaschModel::buildParameterSet(ItemModel*, DimensionModel*) {
}

void RaschModel::setInitialPars(map<Parameter, Matrix<double> >* pair) {
}

void RaschModel::calculateInitialPars() {
}

void RaschModel::successProbability() {
}

const map<Parameter, Matrix<double> *>& RaschModel::getParameterSet() const {
	return this->parameterSet;
}

void RaschModel::setParameterSet(const map<Parameter, Matrix<double> *>& pair) {
	this->parameterSet = parameterSet;
}

RaschModel::~RaschModel() {
	// TODO Auto-generated destructor stub
}

