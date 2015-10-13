/*
 * DimensionModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef DIMENSIONMODEL_H_
#define DIMENSIONMODEL_H_

#include <vector>

using namespace std;

class DimensionModel {
public:
	//Number of dimensions
	int dims = 0;
	// Methods
	virtual int getNumDimensions () = 0;
	virtual vector<int> getDimVector() = 0;
	void setDims(int dim) {
		dims = dim;
	}

	// Destructor
	virtual ~DimensionModel();
};

#endif /* DIMENSIONMODEL_H_ */
