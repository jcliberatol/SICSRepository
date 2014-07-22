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

void RaschModel::successProbability(DimensionModel *) {
}

void RaschModel::setParameterSet(map<Parameter, Matrix<double> *> pair) {
	this->parameterSet = parameterSet;
}

map<Parameter, Matrix<double> *> RaschModel::getParameterSet()  {
	return (this->parameterSet);
}

double RaschModel::getProbability(int node, int item) {
}

RaschModel::~RaschModel() {
	// TODO Auto-generated destructor stub
}

