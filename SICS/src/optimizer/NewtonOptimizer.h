/*
 * NewtonOptimizer.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef NEWTONOPTIMIZER_H_
#define NEWTONOPTIMIZER_H_
#include <optimizer/Optimizer.h>

class NewtonOptimizer : public Optimizer{
public:
	NewtonOptimizer();
	virtual ~NewtonOptimizer();
};

#endif /* NEWTONOPTIMIZER_H_ */