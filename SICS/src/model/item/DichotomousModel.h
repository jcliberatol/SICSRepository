/*
 * DichotomousModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef DICHOTOMOUSMODEL_H_
#define DICHOTOMOUSMODEL_H_

#include <model/item/ItemModel.h>

class DichotomousModel: public ItemModel {
public:
	// Constructor
	DichotomousModel();

	// Methods
	int countCategories();
	  int countItems();

	// Getters and Setters
	PatternMatrix* getDataset();
	void setDataset(PatternMatrix*);
};

#endif /* DICHOTOMOUSMODEL_H_ */
