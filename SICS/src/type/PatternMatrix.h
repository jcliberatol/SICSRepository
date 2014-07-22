/*
 * PatternMatrix.h
 *
 *  Created on: May 30, 2014
 *      Author: mirt
 */

#ifndef PATTERNMATRIX_H_
#define PATTERNMATRIX_H_

#include <boost/dynamic_bitset.hpp>
#include <map>
#include <iostream>
#include <type/DataSet.h>

using namespace std;

class PatternMatrix : public DataSet {

private:
	map<boost::dynamic_bitset<>, long int> matrix;

public:
	//Constructor
	PatternMatrix();

	// Methods
	std::map<boost::dynamic_bitset<>, long int>::const_iterator iterator; //use this when reading in order
	inline void resetIterator(){iterator = matrix.begin();}
	inline bool checkEnd(){return (iterator==matrix.end());}
	inline void iterate(){++iterator;}
	void push(boost::dynamic_bitset<>);//Use this to fill the pattern matrix
	void push(boost::dynamic_bitset<>,int);//Use this to fill the pattern matrix many times with a pattern

	void flush();//Use this to clean matrix

	long int& operator()(boost::dynamic_bitset<>); //Use this to access a specific pattern frecuency and modify it
	friend std::ostream& operator<< (std::ostream &, PatternMatrix &);//Output operator

	//DataSet implementations
	int countItems () const;
	int countIndividuals () const;

	//Destructor
	virtual ~PatternMatrix();
};

#endif /* PATTERNMATRIX_H_ */
