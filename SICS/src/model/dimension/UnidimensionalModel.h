/*
 * UnidimensionalModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef UNIDIMENSIONALMODEL_H_
#define UNIDIMENSIONALMODEL_H_

#include <model/dimension/DimensionModel.h>

class UnidimensionalModel : public DimensionModel {
public:
	// Constructor
	UnidimensionalModel();

	// Methods
	int getNumDimensions ();
	vector<int> getDimVector();
	// Destructor
	virtual ~UnidimensionalModel();
};

#endif /* UNIDIMENSIONALMODEL_H_ */
