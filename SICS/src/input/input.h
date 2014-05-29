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


class input {

private:
	/*
	 * Checks if a file exists, relative to the program
	 */


public:
	input();
	virtual ~input();
	//bool importCSV(char*, GeMatrix<FullStorage<int, RowMajor> >&, int, int);
};

#endif /* INPUT_H_ */
