/*
 * MultiUniDimModel.cpp
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#include <model/dimension/MultiUniDimModel.h>

MultiUniDimModel::MultiUniDimModel() {
	// TODO Auto-generated constructor stub

}

int MultiUniDimModel::getNumDimensions() {
}

vector<double> MultiUniDimModel::getDimVector() {
}

LatentTraitSet* MultiUniDimModel::getLatentTraitSet() const {
	return latentTraitSet;
}

void MultiUniDimModel::setLatentTraitSet(LatentTraitSet* latentTraitSet) {
	this->latentTraitSet = latentTraitSet;
}

MultiUniDimModel::~MultiUniDimModel() {
	// TODO Auto-generated destructor stub
}

