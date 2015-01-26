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
#include <trace/Timer.h>
#include <string>

using namespace std;
using namespace chrono;

/**
 * Trace class for constructing execution and error logs onto files.
 * Create a trace with a file to trace to that file
 * when logging errors onto a trace the error logs onto the file asociated
 * trace object can take ostreams, so any output class of cpp can output to the trace
 * */
class Trace {

	map <string, Timer> timers;
	map <string, string> messages;
	map <string, int> counters;

public:
	const char * filename;

	void startCounter(string id){
		counters[id] = 0;
	}

	void upCount(string id){
		counters[id] = counters[id]+1;
	}

	int readCounter(string id){
		return (counters[id]);
	}

	void storeMessage(string id, string message){
		messages[id]=message;
	}

	string readMessage(string id){
		return (messages[id]);
	}

	void resetTimer(string s){
		timers[s].reset();
	}

	void startTimer(string s){
		timers[s].start();
	}

	void stopTimer(string s){
		timers[s].stop();
	}

	long int timerDuration(string s){
		return (timers[s].totalTime);
	}

	long int dr(string s){
		return (timerDuration(s));
	}
//outputs a matrix
	template<typename T>
	void operator() ( Matrix<T> & message ) {

		ofstream file;
		file.open ( filename, ofstream::app );

		file << message;

		file.close ();
	}
//Outputs a message with a new line
	template<typename T>
	void operator ()(T message) {
		ofstream file;
		file.open ( filename, ofstream::app );

		file << message;

		file << endl;

		file.close ();
	}

//Outputs a message with the specified option
	template<typename T>
	void operator ()(T message, char option) {
		ofstream file;
		file.open ( filename, ofstream::app );

		file << message;

		file.close ();
	}
//Defines the filename for the trace
	Trace( const char * filename ) {
		this->filename = filename;

	}
	~Trace() {
		// TODO Auto-generated destructor stub
	}
	const char* getFilename() const {
		return (filename);
	}

	void setFilename(const char* filename) {
		this->filename = filename;
	}

	void endTrace(){
		ofstream file;
		file.open ( filename, ofstream::app );
		//output of the trace goes here.
		file << endl;

		file.close ();
	}

};
#endif /* TRACE_H_ */
