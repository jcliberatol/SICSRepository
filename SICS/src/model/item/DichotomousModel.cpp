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

int DichotomousModel::countCategories() {
	return (2);
}

DataSet* DichotomousModel::getDataset() {
	return (dataSet);
}
/*
void DichotomousModel::setDataset(DataSet* dataset) {
	this->dataSet = dataset;
}
*/

int DichotomousModel::countItems() {
	PatternMatrix * pM = dynamic_cast<PatternMatrix *>(this->dataSet);
	return (pM->countItems());
}

DichotomousModel::~DichotomousModel() {
	// TODO Auto-generated destructor stub
}

