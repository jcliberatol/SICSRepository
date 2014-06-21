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

class DataSet {
private:
public:
	// Constructor
	DataSet ( );

	// Methods
	vector<int> countItems ();

	// Destructor
	virtual ~DataSet();

};

#endif /* DATASET_H_ */
