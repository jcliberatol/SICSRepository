/*
 * DataSet.h
 *
 *  Created on: 11 Jun 2014
 *      Author: jlgpisa
 */

#ifndef DATASET_H_
#define DATASET_H_

#include <iostream>

/**
 * Skeleton class for the datasets
 * a dataset can contain not only the raw matrices but information about the dataset
 * */
class DataSet
{
	
public:

	// Methods
	virtual   int countItems () const = 0;
	virtual   int countIndividuals () const = 0;

	virtual ~DataSet(){};

};

#endif /* DATASET_H_ */
