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
	ParameterModel *parameterModel = new ThreePLModel();//TODO FIX FOR DECISION OF NEW MODELS
	return (parameterModel);
}

ItemModel* SICSGeneralModel::createItemModel() {
	ItemModel *itemModel = new DichotomousModel();//TODO FIX FOR DECISION OF NEW MODELS
	return (itemModel);
}

DimensionModel* SICSGeneralModel::createDimensionModel() {
	DimensionModel *dimensionModel = new UnidimensionalModel();//TODO FIX FOR DECISION OF NEW MODELS
	return (dimensionModel);
}

SICSGeneralModel::~SICSGeneralModel() {
	// TODO Auto-generated destructor stub
}

