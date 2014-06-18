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
#include <type/Matrix.h>

using namespace std;

class Trace {
	const char * filename;

public:
	Trace( const char * );
	void operator() ( string );
	void operator() ( const char * );
	template<typename T>
	void operator() ( Matrix<T>&);

	virtual ~Trace();
};


template<typename T>
void Trace::operator() ( Matrix<T> & message ) {

	ofstream file;
	file.open ( filename, ofstream::app );

	file << message;

	file.close ();
}


#endif /* TRACE_H_ */
