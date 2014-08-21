/*
 * trace.cpp
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#include "Trace.h"

Trace::Trace( const char * filename ) {

	this->filename = filename;

}


Trace::~Trace() {
	// TODO Auto-generated destructor stub
}

const char* Trace::getFilename() const {
	return filename;
}

void Trace::setFilename(const char* filename) {
	this->filename = filename;
}
