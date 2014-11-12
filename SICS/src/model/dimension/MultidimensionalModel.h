/*
 * MultidimensionalModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef MULTIDIMENSIONALMODEL_H_
#define MULTIDIMENSIONALMODEL_H_

#include <vector>
#include <model/dimension/DimensionModel.h>

using namespace std;

class MultidimensionalModel : public DimensionModel {
public:
	// Constructor
	MultidimensionalModel();

	// Methods
	int getNumDimensions ();
	vector<int> getDimVector();


	// Destructor
	virtual ~MultidimensionalModel();
};

#endif /* MULTIDIMENSIONALMODEL_H_ */
