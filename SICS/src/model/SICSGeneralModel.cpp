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

ParameterModel* SICSGeneralModel::createParameterModel(int model) {

	if (model == Constant::THREE_PL) {
		ParameterModel *parameterModel;
		parameterModel = new ThreePLModel();
		return (parameterModel);
	}

	if (model == Constant::TWO_PL) {
		ParameterModel *parameterModel;
		parameterModel = new TwoPLModel();
		return (parameterModel);
	}

	if (model == Constant::RASCH_A1) {
		ParameterModel *parameterModel;
		// TODO parameterModel = new TwoPLModel();
		return (parameterModel);
	}

	if (model == Constant::RASCH_A_CONSTANT) {
		ParameterModel *parameterModel;
		// TODO parameterModel = new RaschModel();
		return (parameterModel);
	}

	// Default Parameter Model.
	ParameterModel *parameterModel;
	parameterModel = new ThreePLModel();
	return (parameterModel);
}

ItemModel* SICSGeneralModel::createItemModel() {
	ItemModel *itemModel = new DichotomousModel(); //TODO FIX FOR DECISION OF NEW MODELS
	return (itemModel);
}

DimensionModel* SICSGeneralModel::createDimensionModel() {
	DimensionModel *dimensionModel = new UnidimensionalModel(); //TODO FIX FOR DECISION OF NEW MODELS
	return (dimensionModel);
}

SICSGeneralModel::~SICSGeneralModel() {
	// TODO Auto-generated destructor stub
}
