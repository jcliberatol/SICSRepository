/*
 * FisherScoringOptimizer.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef FISHERSCORINGOPTIMIZER_H_
#define FISHERSCORINGOPTIMIZER_H_
#include <core/optimizer/Optimizer.h>

class FisherScoringOptimizer : public Optimizer{
public:
	FisherScoringOptimizer();
	virtual ~FisherScoringOptimizer();
};

#endif /* FISHERSCORINGOPTIMIZER_H_ */
