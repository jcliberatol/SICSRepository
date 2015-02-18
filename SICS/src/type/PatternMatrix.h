/*
 * PatternMatrix.h
 *
 *  Created on: May 30, 2014
 *      Author: mirt
 */

#ifndef PATTERNMATRIX_H_
#define PATTERNMATRIX_H_

#include <map>
#include <iostream>
#include <type/DataSet.h>
#include <vector>

using namespace std;
/**
 * Class for holding binary matrices in the array of patterns of bitsets form.
 * */
class PatternMatrix : public DataSet {

private:
	
public:
	int size;
	map<vector<char>, int> matrix;
	bool** bitset_list;
	int * count_set_bits;
	int * frequency_list;

	//Constructor
	PatternMatrix(int size);

	// Methods
	std::map<vector<char>, int>::const_iterator iterator; /**use this when reading in order*/
	map<vector<char>, int>::const_iterator begin = matrix.begin();
	map<vector<char>, int>::const_iterator end = matrix.end();
	inline void resetIterator(){iterator = matrix.begin();}
	inline bool checkEnd(){return (iterator==matrix.end());}
	inline void iterate(){++iterator;}
	int countBitSet(bool * bitset, int index);
	//inline void values(){cout<<iterator->first<<" values  "<<iterator->second<<endl;}
	inline vector<char> getCurrentBitSet(){return (iterator->first);}
	inline long int getCurrentFrequency(){return (iterator->second);}
	void push(vector<char>);/**Use this to fill the pattern matrix*/
	void push(vector<char>,int);/**Use this to fill the pattern matrix many times with a pattern*/
	int freq(vector<char>);//Frequency
	void flush();/**Use this to clean matrix*/
	void print();


	friend std::ostream& operator<< (std::ostream &, PatternMatrix &);/**Output operator*/

	//DataSet implementations
	int countItems () const;
	int countIndividuals () const;

	bool ** getBitsetList();
	int * getFrequencyList();

	//Destructor
	virtual ~PatternMatrix();
};

#endif /* PATTERNMATRIX_H_ */
