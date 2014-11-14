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

void Model::setModel(ModelFactory * modelFactory, int modelType) {
	parameterModel = modelFactory->createParameterModel(modelType);
	itemModel = modelFactory->createItemModel();
	dimensionModel = modelFactory->createDimensionModel();
}

Model::~Model() {
	delete parameterModel;
	delete itemModel;
	delete dimensionModel;
}
DimensionModel* Model::getDimensionModel()  {
	return (dimensionModel);
}

void Model::setDimensionModel(DimensionModel* dimensionModel) {
	this->dimensionModel = dimensionModel;
}

ItemModel* Model::getItemModel() {
	return (itemModel);
}

void Model::setItemModel(ItemModel* itemModel) {
	this->itemModel = itemModel;
}

ParameterModel* Model::getParameterModel()  {
	return (parameterModel);
}

void Model::setParameterModel(ParameterModel* parameterModel) {
	this->parameterModel = parameterModel;
}

void Model::successProbability(QuadratureNodes* q) {
	parameterModel->successProbability(dimensionModel, q);
}

void Model::buildParameterSet() {
	parameterModel->buildParameterSet(itemModel, dimensionModel);
}
