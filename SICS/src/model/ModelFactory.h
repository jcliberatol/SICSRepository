/*
 * ModelFactory.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef MODELFACTORY_H_
#define MODELFACTORY_H_

class ModelFactory {
public:
	// Constructor
	ModelFactory();

	// Methods
	virtual void createParameterModel() = 0;
	virtual void createItemModel() = 0;
	virtual void createDimensionModel() = 0;

	// Destructor
	virtual ~ModelFactory();
};

#endif /* MODELFACTORY_H_ */
