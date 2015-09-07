/*
 * PolytomousModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef POLYTOMOUSMODEL_H_
#define POLYTOMOUSMODEL_H_

#include <model/item/ItemModel.h>

class PolytomousModel: public ItemModel {
public:
	// Constructor
	PolytomousModel();

	// Methods
	int countCategories();
	unsigned int countItems();

	// Getters and Setters
	PatternMatrix* getDataset() ;
	void setDataset(PatternMatrix* dataset);

	// Destructor
	virtual ~PolytomousModel();
};

#endif /* POLYTOMOUSMODEL_H_ */
