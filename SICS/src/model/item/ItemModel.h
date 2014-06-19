/*
 * ItemModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef ITEMMODEL_H_
#define ITEMMODEL_H_

#include <type/DataSet.h>

class ItemModel {
protected:
	DataSet *dataSet;
public:
	// Methods
	virtual int getCategories () = 0;

	// Getters and Setters
	virtual const DataSet* getDataset() const = 0;
	virtual void setDataset(const DataSet*& dataset) = 0;

	// Destructor
	virtual ~ItemModel();
};

#endif /* ITEMMODEL_H_ */
