/*
 * input.h
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#ifndef INPUT_H_
#define INPUT_H_

#include <boost/filesystem.hpp>
using namespace boost::filesystem;

class input {

private:
	/*
	 * Checks if a file exists, relative to the program
	 */


public:
	input();
	virtual ~input();
	bool checkFile(path file);
};

#endif /* INPUT_H_ */
