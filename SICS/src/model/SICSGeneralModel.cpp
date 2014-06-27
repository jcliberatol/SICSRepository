/*
 * SICSGeneralModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/SICSGeneralModel.h>

SICSGeneralModel::SICSGeneralModel() {
	// TODO Auto-generated constructor stub

}

ParameterModel* SICSGeneralModel::createParameterModel() {
	ParameterModel *parameterModel = new ThreePLModel();
	return parameterModel;
}

ItemModel* SICSGeneralModel::createItemModel() {
	ItemModel *itemModel = new DichotomousModel();
	return itemModel;
}

DimensionModel* SICSGeneralModel::createDimensionModel() {
	DimensionModel *dimensionModel = new UnidimensionalModel();
	return dimensionModel;
}

SICSGeneralModel::~SICSGeneralModel() {
	// TODO Auto-generated destructor stub
}

