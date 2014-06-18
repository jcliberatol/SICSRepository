/*
 * ClassicalEstimation.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef CLASSICALESTIMATION_H_
#define CLASSICALESTIMATION_H_
#include <estimation/Estimation.h>

class ClassicalEstimation : public Estimation {
	Matrix<double> parameters;
public:
	ClassicalEstimation ();
	void estimate ();
	void EM ();

	virtual ~ClassicalEstimation();
};

#endif /* CLASSICALESTIMATION_H_ */
