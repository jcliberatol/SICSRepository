/*
 * QuadratureNodes.cpp
 *
 *  Created on: Oct 16, 2014
 *      Author: jliberato
 */

#include <type/QuadratureNodes.h>

QuadratureNodes::QuadratureNodes() {
	theta = 0;
	weight = 0;
	n=0;
}

int QuadratureNodes::size(){
	return (this->n);
}

QuadratureNodes::QuadratureNodes(Matrix<double>* theta,Matrix<double>* weight) {
	this->theta = theta;
	this->weight = weight;
	this->n=theta->nC();
}

QuadratureNodes::~QuadratureNodes() {
	if (theta) {
		delete theta;
		theta = 0;//TODO remove band aid
	}
	if (weight) {
		delete weight;
		weight = 0;
	}
}

Matrix<double>* QuadratureNodes::getTheta() {
	return (theta);
}

void QuadratureNodes::setTheta(Matrix<double>* theta) {
	this->theta = theta;
	this->n=theta->nC();
}

Matrix<double>* QuadratureNodes::getWeight() {
	return (weight);
}

void QuadratureNodes::setWeight(Matrix<double>* weight) {
	this->weight = weight;
}
