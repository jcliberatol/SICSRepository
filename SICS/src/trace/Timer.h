/*
 * Timer.h
 *
 *  Created on: Dec 23, 2014
 *      Author: jliberato
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <string>

using namespace std;
using namespace chrono;

/*
 * Class Timer
 * Manages timers that can be started an stopped, works in nanoseconds.
 * To Create a timer, reset the timer.
 * Then the timer can be started and stopped by convenience for profiling processes
 * Uses cpp 11 chrono class for measurements
 * Uses the highest resolution clock available
 * Use the totalTime Atribute to measure the counted time.
 */

class Timer {
	typedef high_resolution_clock clock;
	typedef long int millisec;
	typedef time_point<chrono::high_resolution_clock> time;
public:
	millisec totalTime;
	time lastCountPoint;
	bool counting;

	Timer(){
		totalTime = duration_cast<nanoseconds>(nanoseconds(0)).count();
		lastCountPoint = clock::now();
		counting = false;
	}
	void reset(){
		totalTime = duration_cast<nanoseconds>(nanoseconds(0)).count();
		lastCountPoint = clock::now();
		counting = false;
	}

	int start(){
		//In case a timer is already started, starting it will fail
		if (counting == true) {
			lastCountPoint = clock::now();
			return (1);
		}
		else{
			lastCountPoint = clock::now();
			counting = true;
			return (0);
		}

	}

	int stop(){
		if (counting == false) return (1);
		else {
			counting = false;
			totalTime = (totalTime + duration_cast<nanoseconds>(clock::now()-lastCountPoint).count());
			return (0);
		}
	}

	virtual ~Timer(){

	}
};

#endif /* TIMER_H_ */
