/*
 * BayesianEstimation.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef BAYESIANESTIMATION_H_
#define BAYESIANESTIMATION_H_
#include <estimation/Estimation.h>
/**
 * Estimation interface for bayesian methods, currently not implemented
 * */
class BayesianEstimation : public Estimation{
public:
	// Constructor
	BayesianEstimation();

	void estimate ();
	void setModel ( Model * );

	virtual ~BayesianEstimation();
};

#endif /* BAYESIANESTIMATION_H_ */
