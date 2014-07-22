/*
 * Model.h
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <model/parameter/ParameterModel.h>
#include <model/item/ItemModel.h>
#include <model/dimension/DimensionModel.h>
#include <model/ModelFactory.h>

class Model {
	ParameterModel *parameterModel;
	ItemModel *itemModel;
	DimensionModel *dimensionModel;

public:
	// Constructor
	Model();

	// Methods
	void setModel ( ModelFactory * );

	// Getters and Setters
	DimensionModel* getDimensionModel();
	void setDimensionModel(DimensionModel* dimensionModel);

	ItemModel* getItemModel();
	void setItemModel(ItemModel* itemModel);

	ParameterModel* getParameterModel();
	void setParameterModel(ParameterModel* parameterModel);

	// Destructor
	virtual ~Model();
};

#endif /* MODEL_H_ */
