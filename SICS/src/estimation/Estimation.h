/*
 * Estimation.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef ESTIMATION_H_
#define ESTIMATION_H_
#include <model/Model.h>
#include <type/Matrix.h>
#include <optimizer/Optimizer.h>
#include <model/parameter/ThreePLModel.h>
#include <model/parameter/OnePLModel.h>


/**
 * Parent estimation interface for all estimation methods.
 * */
class Estimation {
protected:
	Model *model;

public:
	virtual void estimate () = 0;
	virtual void setModel ( Model * ) = 0;

	// Destructor
	virtual ~Estimation();
};

#endif /* ESTIMATION_H_ */
