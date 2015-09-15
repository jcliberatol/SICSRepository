/*
 * ItemModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef ITEMMODEL_H_
#define ITEMMODEL_H_

#include <type/PatternMatrix.h>

class ItemModel
{

public:

	PatternMatrix *dataSet;
	
	// Methods
	virtual int countCategories() = 0;
	virtual   int countItems() = 0;

	// Getters and Setters
	virtual PatternMatrix* getDataset()  = 0;
	virtual void setDataset(PatternMatrix*); // Dataset must be set before deciding on dimensionality of the model

	// Destructor
	virtual ~ItemModel();
};

#endif /* ITEMMODEL_H_ */
