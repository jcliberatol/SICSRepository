/*
 * EM.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef EM_H_
#define EM_H_
#include <estimation/classical/ClassicalEstimation.h>

class EMEstimation : public ClassicalEstimation{
public:
	EMEstimation();
	virtual ~EMEstimation();
};

#endif /* EM_H_ */
