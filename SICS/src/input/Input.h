/*
 * input.h
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#ifndef INPUT_H_
#define INPUT_H_

#include <map>
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include <iostream>
#include <ctype.h>
#include <fstream>
#include <type/PatternMatrix.h>
#include <trace/Trace.h>

using namespace std;


class Input {

private:
	/*
	 * Checks if a file exists, relative to the program
	 */
	char del = ',';

public:
	Input();
	virtual ~Input();
	bool importCSV ( char*, PatternMatrix&, unsigned int, unsigned int );

	char getDel() const;
	void setDel(char);

};

#endif /* INPUT_H_ */
