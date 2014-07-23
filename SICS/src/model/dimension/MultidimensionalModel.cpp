/*
 * MultidimensionalModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include "MultidimensionalModel.h"

MultidimensionalModel::MultidimensionalModel() {
	// TODO Auto-generated constructor stub

}

vector<int> MultidimensionalModel::getDimVector() {
}

MultidimensionalModel::~MultidimensionalModel() {
	// TODO Auto-generated destructor stub
}

LatentTraitSet* MultidimensionalModel::getLatentTraitSet() const {
	return latentTraitSet;
}

void MultidimensionalModel::setLatentTraitSet(
		LatentTraitSet* latentTraitSet) {
	this->latentTraitSet = latentTraitSet;
}

int MultidimensionalModel::getNumDimensions() {
}

