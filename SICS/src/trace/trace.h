/*
 * trace.h
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */
#include <iostream>
#ifndef TRACE_H_
#define TRACE_H_

class trace {
private:
public:
	trace();
	void log(std::ostream &);
	virtual ~trace();
};

#endif /* TRACE_H_ */
