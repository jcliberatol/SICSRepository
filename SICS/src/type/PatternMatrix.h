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
/**
 * Class for holding binary matrices in the array of patterns of bitsets form.
 * */
class PatternMatrix : public DataSet {

private:
	
public:
	int size;
	//map<boost::dynamic_bitset<>,int> matrix;
	map<bool*, int> matrix;
	//Constructor
	PatternMatrix(int size);

	// Methods
	//std::map<boost::dynamic_bitset<>, int>::const_iterator iterator; /**use this when reading in order*/
	std::map<bool*, int>::const_iterator iterator; /**use this when reading in order*/
	inline void resetIterator(){iterator = matrix.begin();}
	inline bool checkEnd(){return (iterator==matrix.end());}
	inline void iterate(){++iterator;}
	inline void values(){cout<<iterator->first<<" values  "<<iterator->second<<endl;}
	//inline boost::dynamic_bitset<> getCurrentBitSet(){return (iterator->first);}
	inline bool* getCurrentBitSet(){return (iterator->first);}
	inline long int getCurrentFrequency(){return (iterator->second);}
//	void push(boost::dynamic_bitset<>);/**Use this to fill the pattern matrix*/
//	void push(boost::dynamic_bitset<>,int);/**Use this to fill the pattern matrix many times with a pattern*/
	void push(bool*);/**Use this to fill the pattern matrix*/
	void push(bool*,int);/**Use this to fill the pattern matrix many times with a pattern*/

	void flush();/**Use this to clean matrix*/

	//int& operator()(boost::dynamic_bitset<>); /**Use this to access a specific pattern frecuency and modify it*/
	friend std::ostream& operator<< (std::ostream &, PatternMatrix &);/**Output operator*/

	//DataSet implementations
	int countItems () const;
	int countIndividuals () const;

	//Destructor
	virtual ~PatternMatrix();
};

#endif /* PATTERNMATRIX_H_ */
