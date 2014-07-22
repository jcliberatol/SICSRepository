/*
 * LatentTraitSet.h
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#ifndef LATENTTRAITSET_H_
#define LATENTTRAITSET_H_

#include <type/Matrix.h>

class LatentTraitSet {
	Matrix<double> *theta;
	Matrix<double> *weight;
public:
	// Constructor
	LatentTraitSet();
	
	// Getters and Setters  
	Matrix<double>* getTheta() const;
	void setTheta(Matrix<double>* theta);
	Matrix<double>* getWeight() const;
	void setWeight(Matrix<double>* weight);

	// Destructor
	virtual ~LatentTraitSet();

};

#endif /* LATENTTRAITSET_H_ */
