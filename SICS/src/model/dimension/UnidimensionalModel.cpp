/*
 * UnidimensionalModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/dimension/UnidimensionalModel.h>

UnidimensionalModel::UnidimensionalModel() {

}

int UnidimensionalModel::getNumDimensions() {
	return (1);
}

vector<int> UnidimensionalModel::getDimVector() {
	vector<int> dimVect;

	int dimension = 0;//latentTraitSet->getTheta()->nC();
	dimVect.push_back(dimension);

	return (dimVect);
}

UnidimensionalModel::~UnidimensionalModel(){

}

