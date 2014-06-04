/*
 * trace.h
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#ifndef TRACE_H_
#define TRACE_H_

#include <string>
#include <fstream>

using namespace std;

class Trace {
	const char * filename;

public:
	Trace( const char * );
	void operator() ( string );
	void operator() ( const char * );
	void operator() ( ostream & );

	virtual ~Trace();
};

#endif /* TRACE_H_ */
