/*
 * DichotomousModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef DICHOTOMOUSMODEL_H_
#define DICHOTOMOUSMODEL_H_

#include <model/item/ItemModel.h>
#include <type/PatternMatrix.h>

class DichotomousModel: public ItemModel {
public:
	// Constructor
	DichotomousModel();

	// Methods
	int countCategories();
	int countItems();

	// Getters and Setters
	DataSet* getDataset();

	// Destructor
	virtual ~DichotomousModel();
};

#endif /* DICHOTOMOUSMODEL_H_ */
