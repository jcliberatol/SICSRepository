/*
 * BFGSOptimizer.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef BFGSOPTIMIZER_H_
#define BFGSOPTIMIZER_H_
#include <optimizer/Optimizer.h>

class BFGSOptimizer : public Optimizer{
public:
	BFGSOptimizer();
	virtual ~BFGSOptimizer();
};

#endif /* BFGSOPTIMIZER_H_ */
