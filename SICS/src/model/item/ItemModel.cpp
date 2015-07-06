/*
 * ItemModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/item/ItemModel.h>
#include <iostream>

ItemModel::~ItemModel() {}

void ItemModel::setDataset(PatternMatrix* dataset) { this->dataSet = dataset; }
