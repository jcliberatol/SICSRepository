/*
 * EM.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef EM_H_
#define EM_H_
#include <estimation/classical/ClassicalEstimation.h>

class EM : public ClassicalEstimation{
public:
	EM();
	virtual ~EM();
};

#endif /* EM_H_ */
