/*
 * Model.cpp
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#include "Model.h"

Model::Model ( ) {
	// TODO Auto-generated constructor stub

}

void Model::setModel ( ModelFactory * modelFactory ) {
}

Model::~Model() {
	// TODO Auto-generated destructor stub
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

const ParameterModel* Model::getParameterModel() const {
	return parameterModel;
}

void Model::setParameterModel(ParameterModel* parameterModel) {
	this->parameterModel = parameterModel;
}
