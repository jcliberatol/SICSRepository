/*
 * ItemModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/item/ItemModel.h>
#include <iostream>

ItemModel::~ItemModel()
{
	// TODO
	// if(dataSet != NULL)
	// 	delete dataSet;
}

void ItemModel::setDataset(DataSet* dataset) { this->dataSet = dataset; }
