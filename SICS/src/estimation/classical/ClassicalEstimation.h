/*
 * ClassicalEstimation.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef CLASSICALESTIMATION_H_
#define CLASSICALESTIMATION_H_
#include <estimation/Estimation.h>

/**
 * Classical estimation class all the classical estimation methods must derive from this class
 * this class is an interface.
 * */
class ClassicalEstimation : public Estimation {
public:
	ClassicalEstimation ();

	virtual void estimate () = 0;/** Estimation method for the item parameters , after executed, the item parameters must be equal to the best parameters for the item according to the model*/
	virtual void setModel ( Model * ) = 0;/** Sets the model for the classical estimation to use.*/

	virtual ~ClassicalEstimation();
};

#endif /* CLASSICALESTIMATION_H_ */
