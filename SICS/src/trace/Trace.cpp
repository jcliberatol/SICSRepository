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

void Trace::operator() ( string message ) {

	ofstream file;
	file.open ( filename, ofstream::app );

	file << message << endl;

	file.close ();

}

void Trace::operator() ( const char * message ) {

	ofstream file;
	file.open ( filename, ofstream::app );

	file << message << endl;

	file.close ();
}


Trace::~Trace() {
	// TODO Auto-generated destructor stub
}
