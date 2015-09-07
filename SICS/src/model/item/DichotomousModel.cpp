/*
 * DichotomousModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/item/DichotomousModel.h>

DichotomousModel::DichotomousModel() {}

int DichotomousModel::countCategories() { return (2); }

PatternMatrix* DichotomousModel::getDataset() { return (dataSet); }

void DichotomousModel::setDataset(PatternMatrix* dataset) { this->dataSet = dataset; }


unsigned int DichotomousModel::countItems()
{
	unsigned int result = this->dataSet->countItems();
	return (result);
}

