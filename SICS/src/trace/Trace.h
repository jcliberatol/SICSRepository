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
#include <chrono>

using namespace std;

/**
 * Trace class for constructing execution and error logs onto files.
 * Create a trace with a file to trace to that file
 * when logging errors onto a trace the error logs onto the file asociated
 * trace object can take ostreams, so any output class of cpp can output to the trace
 * */
class Trace {
	const char * filename;
	chrono::time_point<std::chrono::system_clock> start, end;
	chrono::duration<double> timeStepE;
	chrono::duration<double> timeStepM;
	chrono::duration<double> timeSetInitialValues;



public:

	template<typename T>
	void operator() ( Matrix<T> & message ) {

		ofstream file;
		file.open ( filename, ofstream::app );

		file << message;

		file.close ();
	}

	template<typename T>
	void operator ()(T message) {
		ofstream file;
		file.open ( filename, ofstream::app );

		file << message;

		file << endl;

		file.close ();
	}


	Trace( const char * filename ) {

		this->filename = filename;

	}
	~Trace() {
		// TODO Auto-generated destructor stub
	}

	const char* getFilename() const {
		return filename;
	}

	void setFilename(const char* filename) {
		this->filename = filename;
	}

	void startTimingMeasure(){
		this->start = chrono::system_clock::now();
	}
	/**
	 *  Type
	 *  	0 = stepE
	 *  	1 = stepM
	 *  	2 = setInitialValues
	 */
	void finishTimingMeasure(int type){
		this->end = chrono::system_clock::now();
		switch (type) {
		case 0:
			this->timeStepE += end - start;
			break;
		case 1:
			this-> timeStepM += end - start;
			break;
		case 2:
			this ->timeSetInitialValues += end -start;
			break;
		default:
			break;
		}
	}

	void endTrace(){
		ofstream file;
		file.open ( filename, ofstream::app );
		chrono::duration<double> total= timeSetInitialValues + timeStepE + timeStepM;
		file << "Step E: "<<timeStepE.count() <<"s "<< (timeStepE/total)*100 << "%\n";
		file << "Step M: "<<timeStepM.count() <<"s "<< (timeStepM/total)*100 << "%\n";
		file << "Set initial values: "<<timeSetInitialValues.count() <<"s "<< (timeSetInitialValues/total)*100 << "%\n";

		file << endl;

		file.close ();
	}

};
#endif /* TRACE_H_ */
