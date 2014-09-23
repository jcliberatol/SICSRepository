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
#include <type/Matrix.h>
#include <trace/Trace.h>

using namespace std;

/**
 * Class that is in charge of taking OS files, streams and other sources for inputting data into the software suite
 * */
class Input {

private:
	/*
	 * Checks if a file exists, relative to the program
	 */
	char del;

public:
	Input();
	virtual ~Input();
	bool importCSV ( char*, PatternMatrix&, unsigned int, unsigned int );/**Imports binary matrices from a csv*/
	bool importCSV ( char*, Matrix<double>&, unsigned int , unsigned int);/**Imports generic type matrices from a csv*/
	char getDel() const;/**Gets the delimitier used for inputting*/
	void setDel(char);/**Sets the delimitier for inputting text matrices*/

};

#endif /* INPUT_H_ */
