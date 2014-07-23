/*
 * MultidimensionalModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef MULTIDIMENSIONALMODEL_H_
#define MULTIDIMENSIONALMODEL_H_

#include <vector>
#include <type/LatentTraitSet.h>
#include <model/dimension/DimensionModel.h>

using namespace std;

class MultidimensionalModel : public DimensionModel {
public:
	// Constructor
	MultidimensionalModel();

	// Methods
	int getNumDimensions ();
	vector<int> getDimVector();

	// Getters and Setters
	LatentTraitSet* getLatentTraitSet() const;
	void setLatentTraitSet(LatentTraitSet* latentTraitSet);

	// Destructor
	virtual ~MultidimensionalModel();
};

#endif /* MULTIDIMENSIONALMODEL_H_ */
