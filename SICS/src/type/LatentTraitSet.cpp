/*
 * LatentTraitSet.cpp
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#include <type/LatentTraitSet.h>

LatentTraitSet::LatentTraitSet() {
	// TODO Auto-generated constructor stub

}

LatentTraitSet::~LatentTraitSet() {
	// TODO Auto-generated destructor stub
}

Matrix<double>* LatentTraitSet::getTheta() const {
	return theta;
}

void LatentTraitSet::setTheta(Matrix<double>* theta) {
	this->theta = theta;
}

Matrix<double>* LatentTraitSet::getWeight() const {
	return weight;
}

void LatentTraitSet::setWeight(Matrix<double>* weight) {
	this->weight = weight;
}
