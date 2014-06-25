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

	if ( typeid(*itemModel) == typeid(DichotomousModel) ) {

		if ( typeid(*dimensionModel) == typeid (UnidimensionalModel) ) {

			int items = itemModel->countItems();

			parameterSet[a] = new Matrix<double> (1,items);
			parameterSet[d] = new Matrix<double> (1,items);
			parameterSet[c] = new Matrix<double> (1,items);

		}

		else if ( typeid(*dimensionModel) == typeid (MultidimensionalModel) ) {
			// TODO: Dichotomous Multidimensional
		}

		else if ( typeid(*dimensionModel) == typeid (MultiUniDimModel) ) {
			// TODO: Dichotomous MultiUniDimensional
		}

	}

	else if ( typeid(*dimensionModel) == typeid(PolytomousModel) ) {
		// TODO: Polytomous Model for Unidimensional, Multidimensional and MultiUni
	}

}

void ThreePLModel::setInitialPars(
	map<Parameter, Matrix<double> >* parameterSet) {
}

void ThreePLModel::calculateInitialPars() {
}

void ThreePLModel::successProbability() {
}

const map<Parameter, Matrix<double> *>& ThreePLModel::getParameterSet() const {
return this->parameterSet;
}

void ThreePLModel::setParameterSet(const map<Parameter, Matrix<double> *>& parameterSet) {
	this->parameterSet = parameterSet;
}

ThreePLModel::~ThreePLModel() {
// TODO Auto-generated destructor stub
}

