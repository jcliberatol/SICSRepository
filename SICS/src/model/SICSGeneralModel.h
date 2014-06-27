/*
 * SICSGeneralModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef SICSGENERALMODEL_H_
#define SICSGENERALMODEL_H_

#include <model/ModelFactory.h>
#include <model/parameter/ThreePLModel.h>
#include <model/item/DichotomousModel.h>
#include <model/dimension/UnidimensionalModel.h>

class SICSGeneralModel : public ModelFactory {
public:
	// Constructor
	SICSGeneralModel();

	// Methods
	ParameterModel *createParameterModel();
	ItemModel *createItemModel();
	DimensionModel *createDimensionModel();

	// Destructor
	virtual ~SICSGeneralModel();
};

#endif /* SICSGENERALMODEL_H_ */
