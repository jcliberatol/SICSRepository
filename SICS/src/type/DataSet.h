/*
 * DataSet.h
 *
 *  Created on: 11 Jun 2014
 *      Author: jlgpisa
 */

#ifndef DATASET_H_
#define DATASET_H_

#include <type/PatternMatrix.h>
#include <vector>

using namespace std;

enum DatasetType { PATTERN_MATRIX };

class DataSet {
private:
	DatasetType datasetType;;
	PatternMatrix *dataSet;
public:
	DataSet ( PatternMatrix* );
	DatasetType getDatasetType() const;
	void setDatasetType(DatasetType datasetType);
	vector<int> getDimensions();
	//Dataset must define methods that allow the initial values algorithms to easily access to parameters they may need about the dataset,
	//Size of the dataset, total sum of the patterns , sum of a specific column, accessors to specific columns or rows, etc.
	virtual ~DataSet();

};

#endif /* DATASET_H_ */
