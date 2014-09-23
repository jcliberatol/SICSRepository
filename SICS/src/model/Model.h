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
/**
 * Model class that holds the structures for the IRT models
 * can vary across parameters, items and dimensions
 * includes suport for dichotomic and polytomic models
 * multidimensional and singledimensional models
 * future suport for multiscale and longitudinal models can be implemented.
 * */
class Model {
	ParameterModel *parameterModel;
	ItemModel *itemModel;
	DimensionModel *dimensionModel;

public:
	// Constructor
	Model();

	// Methods
	void setModel ( ModelFactory * );
	void successProbability ();
	void buildParameterSet ();

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
