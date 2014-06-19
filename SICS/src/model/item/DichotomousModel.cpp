/*
 * DichotomousModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/item/DichotomousModel.h>

DichotomousModel::DichotomousModel() {
	// TODO Auto-generated constructor stub

}

int DichotomousModel::getCategories() {
	return 2;
}

const DataSet* DichotomousModel::getDataset() const {
	return dataSet;
}

void DichotomousModel::setDataset(DataSet* dataset) {
	this->dataSet = dataSet;
}

DichotomousModel::~DichotomousModel() {
	// TODO Auto-generated destructor stub
}

