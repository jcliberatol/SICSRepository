/*
 * UnidimensionalModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/dimension/UnidimensionalModel.h>

UnidimensionalModel::UnidimensionalModel() {
	latentTraitSet = new LatentTraitSet();
}

int UnidimensionalModel::getNumDimensions() {
	return (1);
}

vector<int> UnidimensionalModel::getDimVector() {
	vector<int> dimVect;

	int dimension = latentTraitSet->getTheta()->nC();
	dimVect.push_back(dimension);

	return (dimVect);
}

LatentTraitSet* UnidimensionalModel::getLatentTraitSet() const {
	return (latentTraitSet);
}

void UnidimensionalModel::setLatentTraitSet(LatentTraitSet* latentTraitSet) {
	this->latentTraitSet = latentTraitSet;
}

UnidimensionalModel::~UnidimensionalModel() {
	if (latentTraitSet!=NULL) {
		delete latentTraitSet;
	}
}

