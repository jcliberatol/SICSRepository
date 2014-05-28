/*
 * input.cpp
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#include "input.h"

input::input() {
	// TODO Auto-generated constructor stub

}

input::~input() {
	// TODO Auto-generated destructor stub
}

bool input::checkFile(path file){
	return exists(file);
}
