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

const DimensionModel& Model::getDimensionModel() const {
	return dimensionModel;
}

void Model::setDimensionModel(const DimensionModel& dimensionModel) {
	this->dimensionModel = dimensionModel;
}

const ItemModel& Model::getItemModel() const {
	return itemModel;
}

void Model::setItemModel(const ItemModel& itemModel) {
	this->itemModel = itemModel;
}

const ParameterModel& Model::getParameterModel() const {
	return parameterModel;
}

void Model::setParameterModel(const ParameterModel& parameterModel) {
	this->parameterModel = parameterModel;
}

Model::~Model() {
	// TODO Auto-generated destructor stub
}
