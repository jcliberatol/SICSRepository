/*
 * LatentTraitEstimation.h
 *
 *  Created on: Oct 16, 2014
 *      Author: jliberato
 */

#ifndef LATENTTRAITESTIMATION_H_
#define LATENTTRAITESTIMATION_H_
#include <model/Model.h>

class LatentTraitEstimation {
public:
	LatentTraitEstimation();
	virtual ~LatentTraitEstimation();
	void setModel(Model* model);

private:
	Model* model;
};

#endif /* LATENTTRAITESTIMATION_H_ */
