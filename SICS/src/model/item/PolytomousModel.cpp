/*
 * PolytomousModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/item/PolytomousModel.h>

PolytomousModel::PolytomousModel() {
	// TODO Auto-generated constructor stub

}

int PolytomousModel::countCategories() {
	return (0);
}

PatternMatrix* PolytomousModel::getDataset()  {
	return (dataSet);
}

void PolytomousModel::setDataset(PatternMatrix* dataset) {
	this->dataSet = dataset;
}

int PolytomousModel::countItems() {
	return (0);
}

PolytomousModel::~PolytomousModel() {
	// TODO Auto-generated destructor stub
}

