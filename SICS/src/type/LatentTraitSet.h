/*
 * LatentTraitSet.h
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#ifndef LATENTTRAITSET_H_
#define LATENTTRAITSET_H_

#include <type/Matrix.h>

/*
 * Notice this class models the fictious latent trait sets used to estimate the item parameters
 * and not the actual individual latent trait sets
 */
class LatentTraitSet {
	Matrix<double> *theta;
	Matrix<double> *weight;
public:
	// Constructor
	LatentTraitSet();
	
	// Getters and Setters  
	Matrix<double>* getTheta();
	void setTheta(Matrix<double>* theta);
	Matrix<double>* getWeight();
	void setWeight(Matrix<double>* weight);

	// Destructor
	virtual ~LatentTraitSet();

};

#endif /* LATENTTRAITSET_H_ */
