/*
 * DataSet.cpp
 *
 *  Created on: 11 Jun 2014
 *      Author: jlgpisa
 */

#include <type/DataSet.h>

DataSet::DataSet( PatternMatrix *pM ) {

	dataSet = pM;
	datasetType = PATTERN_MATRIX;

}

DatasetType DataSet::getDatasetType() const {
	return datasetType;
}

void DataSet::setDatasetType(DatasetType datasetType) {
	this->datasetType = datasetType;
}

std::vector<int> DataSet::getDimensions() {
	std::vector<int> dimensions;

	//Number of items
	dimensions.push_back(dataSet->getBitsetLength());
	//TODO: Number of individuals: dimensions.push_back();

	return dimensions;
}

DataSet::~DataSet() {
	// TODO Auto-generated destructor stub
}

