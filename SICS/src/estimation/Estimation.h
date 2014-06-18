/*
 * Estimation.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef ESTIMATION_H_
#define ESTIMATION_H_
#include <type/DataSet.h>
#include <model/ModelFactory.h>
#include <model/Model.h>
#include <type/Matrix.h>

class Estimation {
private:
	Model *model;

public:
	// Constructor
	Estimation();

	virtual void setModel ( Model * ) = 0;
	virtual void estimate () = 0;

	// Destructor
	virtual ~Estimation();
};

#endif /* ESTIMATION_H_ */
