/*
 * LatentTraitSet.cpp
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#include <type/LatentTraitSet.h>

LatentTraitSet::LatentTraitSet() {
	theta = NULL;
	weight = NULL;
}

Matrix<double>* LatentTraitSet::getTheta() {
	return (theta);
}

void LatentTraitSet::setTheta(Matrix<double>* theta) {
	this->theta = theta;
}

Matrix<double>* LatentTraitSet::getWeight() {
	return (weight);
}

void LatentTraitSet::setWeight(Matrix<double>* weight) {
	this->weight = weight;
}

LatentTraitSet::~LatentTraitSet() {
	if (theta) {
		delete theta;
		theta = 0;
	}
	if (weight) {
		delete weight;
		weight = 0;
	}
}
