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
map<Parameter, Matrix<double> *> TwoPLModel::getParameterSet()  {
	return (this->parameterSet);
}

double TwoPLModel::successProbability(double theta, double a, double b, double d) {
	long double exponential = (Constant::D_CONST*((a*theta)-(a*b)));
	if ( exponential > Constant::MAX_EXP ) {
		exponential = Constant::MAX_EXP;
	}
	
	else if ( exponential < -(Constant::MAX_EXP*1.0) ) {
		exponential = -Constant::MAX_EXP;
	}
	
	exponential = exp(-exponential);
	
	return ( 1 / ( 1 + exponential ) );
}

void TwoPLModel::setParameterSet(map<Parameter, Matrix<double> *> pair) {
	this->parameterSet = parameterSet;
}

double TwoPLModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}

map<Parameter, Matrix<double> *> TwoPLModel::getParameterSet() {
	return (this->parameterSet);
}

TwoPLModel::~TwoPLModel() {
	// TODO Auto-generated destructor stub
}

