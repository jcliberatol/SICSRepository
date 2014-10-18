/*
 * QuadratureNodes.h
 *
 *  Created on: Oct 16, 2014
 *      Author: jliberato
 */

#ifndef QUADRATURENODES_H_
#define QUADRATURENODES_H_
#include <type/Matrix.h>

class QuadratureNodes {
private:
	Matrix<double> *theta;
	Matrix<double> *weight;
	int n;
public:
	// Getters and Setters
	Matrix<double>* getTheta();
	void setTheta(Matrix<double>* theta);
	Matrix<double>* getWeight();
	void setWeight(Matrix<double>* weight);
	QuadratureNodes();
	QuadratureNodes(Matrix<double>*,Matrix<double>*);
	int size();
	virtual ~QuadratureNodes();
};

#endif /* QUADRATURENODES_H_ */
