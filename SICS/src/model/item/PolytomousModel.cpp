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
}

 DataSet* PolytomousModel::getDataset()  {
	return dataSet;
}

void PolytomousModel::setDataset(DataSet* dataset) {
	this->dataSet = dataset;
}

int PolytomousModel::countItems() {
}

PolytomousModel::~PolytomousModel() {
	// TODO Auto-generated destructor stub
}

