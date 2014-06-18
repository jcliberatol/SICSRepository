/*
 * PatternMatrix.h
 *
 *  Created on: May 30, 2014
 *      Author: mirt
 */
#include <boost/dynamic_bitset.hpp>
#include <map>
#include <iostream>



#ifndef PATTERNMATRIX_H_
#define PATTERNMATRIX_H_

class PatternMatrix {

private:
	std::map<boost::dynamic_bitset<>, long int> matrix;

public:
	PatternMatrix();
	std::map<boost::dynamic_bitset<>, long int>::const_iterator iterator; //use this when reading in order

	void push(boost::dynamic_bitset<>);//Use this to fill the pattern matrix
	void push(boost::dynamic_bitset<>,int);//Use this to fill the pattern matrix many times with a pattern

	void flush();//Use this to clean matrix

	long int& operator()(boost::dynamic_bitset<>); //Use this to access a specific pattern frecuency and modify it
	friend std::ostream& operator<< (std::ostream &, PatternMatrix &);//Output operator

	int getBitsetLength() const;

	virtual ~PatternMatrix();
};

#endif /* PATTERNMATRIX_H_ */