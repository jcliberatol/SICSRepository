/*
 * TwoPLModel.cpp
 *
 *  Created on: 18 Jun 2014
 *  Updated on: 22 Oct 2014
 *      Author: cesandvalp
 */

#include <model/parameter/TwoPLModel.h>

TwoPLModel::TwoPLModel() {
	// TODO Auto-generated constructor stub

}

void TwoPLModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {
			if (typeid(*itemModel) == typeid(DichotomousModel)) {
				if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {
					int items = itemModel->countItems();
					parameterSet[a] = new Matrix<double>(1, items);
					parameterSet[b] = new Matrix<double>(1, items);
					parameterSet[d] = new Matrix<double>(1, items);
				}
				
				else if (typeid(*dimensionModel) == typeid(MultidimensionalModel)) {
					// TODO: Dichotomous Multidimensional
				}
				
				else if (typeid(*dimensionModel) == typeid(MultiUniDimModel)) {
					// TODO: Dichotomous MultiUniDimensional
				}
			}
			
			else if (typeid(*dimensionModel) == typeid(PolytomousModel)) {
				// TODO: Polytomous Model for Unidimensional, Multidimensional and MultiUni
			}
}

void TwoPLModel::successProbability(DimensionModel *dimensionModel,  QuadratureNodes *quadNodes) {
	int q = 0;
	double a_d, b_d, d_d, theta_d; // d stands from "double"
	
	if (dimensionModel != NULL) {
		q = quadNodes->size();
	}
	
	if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {
		int It = parameterSet[a]->nC();
		if (probabilityMatrix == NULL) {
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q, It);
		}
		for (int k = 0; k < q; k++) {
			for (int i = 0; i < It; i++) {
				// 2PL Success Probability Function
				theta_d = (*quadNodes->getTheta())(0, k);
				a_d = (*parameterSet[a])(0, i);
				b_d = (*parameterSet[b])(0, i);
				d_d = (*parameterSet[d])(0, i);
				double p_d = successProbability(theta_d, a_d, b_d, d_d);
				(*probabilityMatrix)(k, i) = p_d;
			}
		}
	}
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

