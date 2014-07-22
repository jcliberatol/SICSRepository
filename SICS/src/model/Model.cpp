/*
 * Model.cpp
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#include "Model.h"

Model::Model ( ) {
	parameterModel = NULL;
	itemModel = NULL;
	dimensionModel = NULL;
}

void Model::setModel(ModelFactory * modelFactory) {
	parameterModel = modelFactory->createParameterModel();
	itemModel = modelFactory->createItemModel();
	dimensionModel = modelFactory->createDimensionModel();
}

Model::~Model() {
	delete parameterModel;
	delete itemModel;
	delete dimensionModel;
}

const DimensionModel* Model::getDimensionModel() const {
	return dimensionModel;
}

void Model::setDimensionModel(DimensionModel* dimensionModel) {
	this->dimensionModel = dimensionModel;
}

const ItemModel* Model::getItemModel() const {
	return itemModel;
}

void Model::setItemModel(ItemModel* itemModel) {
	this->itemModel = itemModel;
}

ParameterModel* Model::getParameterModel()  {
	return parameterModel;
}

void Model::setParameterModel(ParameterModel* parameterModel) {
	this->parameterModel = parameterModel;
}
