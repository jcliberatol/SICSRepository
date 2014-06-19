/*
 * UnidimensionalModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/dimension/UnidimensionalModel.h>

UnidimensionalModel::UnidimensionalModel() {
	// TODO Auto-generated constructor stub

}

int UnidimensionalModel::getNumDimensions() {
}

vector<double> UnidimensionalModel::getDimVector() {
}

const LatentTraitSet* UnidimensionalModel::getLatentTraitSet() const {
	return latentTraitSet;
}

void UnidimensionalModel::setLatentTraitSet(LatentTraitSet* latentTraitSet) {
	this->latentTraitSet = latentTraitSet;
}

UnidimensionalModel::~UnidimensionalModel() {
	// TODO Auto-generated destructor stub
}

