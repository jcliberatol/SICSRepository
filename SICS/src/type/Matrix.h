/*
 * Matrix.h
 *
 *  Created on: May 28, 2014
 *      Author: mirt
 */
#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <climits>

using namespace std;


/**
 * Supports the Matrix Structures, with indexing on any type,
 *  has special operations for general matrix use, is the class used by all the BLAS interface methods
 * , please use this class in the package whenever posible.
 * the class is a templated class supporting any type for components of the matrix, the support of this class is
 * made for two dimension matrices and one dimension vectors, overloaded operators are the parenthesis for one or two
 * dimensional indexing, and the output operator, the class also supports transposed and symmetric matrices
 * and when declaring matrices they can be declared either empty, reseted to a value, as an identity, or as a
 * random matrix. Some very fast methods are implemented for specific matrix sizes. Use this with care.
 */
template<typename T>
class Matrix;

/**
 * Overloads the << operator for printing two dimensional matrices, use for debugging or outputting to files
 * returns an ostream object
 */
template<typename T>
ostream& operator<<(ostream &, Matrix<T> &);

template<class T>
class Matrix {
private:
	unsigned int nCol;
	unsigned int nRow;
	T m(char);
	T get3x3determinant ();

public:
	//bool transposed;
	bool symmetric;
	T **memory;
	unsigned int ld;
	static char del;
	Matrix(T**, unsigned int, unsigned int);/*Importing a 2D array*/
	Matrix(); /**Empty object*/
	Matrix(unsigned int, unsigned int); /**Two dimensional Matrix Constructor allocates memory*/
	Matrix(Matrix<T>&); /**Copy constructor*/
	Matrix(char I, unsigned int size); /**Create special kinds of matrices (dense identity)*/
	void reset();/**Reset method, puts all entries in zeros*/
	void transpose ();/**Transposes the matrix, notice it does not perform memory transpose, only index transpose*/
	void copy(Matrix<T>&);/**Copy constructor*/
	unsigned int nR(); /** Returns number of rows */
	unsigned int nC(); /** Returns number of columns */
	void setIndex(unsigned int, unsigned int, T);
	T getIndex(unsigned int, unsigned int);
	T sum(); /** Returns the sum of all objects */
	inline T & operator()(const unsigned int nCol, const unsigned int nRow); /** Accessing operator for a element */
	T & operator()(const unsigned int element); /**Accessing operator for a element */
	friend ostream& operator<<<T>(ostream &, Matrix<T> &); /** Output operator */
	bool isSymmetric() const; /** Symmetry flag for optimizations */
	void setSymmetric(bool symmetric); /** Set to true the symmetry flag */
	virtual ~Matrix();
};

template<class T>
char Matrix<T>::del = ' ';

template<class T>
T Matrix<T>::sum()
{
  T sum = 0;

  for (unsigned int i = 0; i < nRow; i++)
    for (unsigned int j = 0; j < nCol; j++)
      sum += memory[i];

  return (sum);
}

template<class T>
unsigned int Matrix<T>::nR() { return (nRow); }

template<class T>
unsigned int Matrix<T>::nC() { return (nCol); }

template<class T>
inline void Matrix<T>::setIndex(unsigned int r, unsigned int c, T value)
{ memory[r][c] = value; }

template<class T>
inline T Matrix<T>::getIndex(unsigned int r, unsigned int c)
{ return memory[r][c]; }

template<class T>
Matrix<T>::Matrix() {
	nCol = 0;
	nRow = 0;
	ld = 0;
	memory = NULL;
	//transposed = false;
	symmetric = false;
}

template<class T>
Matrix<T>::Matrix(Matrix<T>& a) {
	memory = NULL;
	copy(a);
}

template<class T>
Matrix<T>::Matrix(unsigned int r, unsigned int c) {
	nCol = c;
	nRow = r;
	//transposed = false;
	memory = new T*[r];
	for(unsigned int i = 0; i < r; i++)
	  memory[i] = new T[c];
	
	symmetric = false;
	ld = c;
}

template<class T>
Matrix<T>::Matrix(T** mem , unsigned int r, unsigned int c) {
	nCol = c;
	nRow = r;
	//transposed = false;
	memory = new T*[r];
	for(unsigned int i = 0; i < r; i++)
	  memory[i] = new T[c];
	symmetric = false;
	ld = c;

	for (unsigned int var = 0; var < r; ++var) {
		for (unsigned int j = 0; j < c; ++j) {
			memory[var][j] = mem[var][j];
		}
	}
}

template<class T>
inline T & Matrix<T>::operator()(const unsigned int r, const unsigned int c) {
		return (memory[r][c]);

}

template<class T>
inline bool Matrix<T>::isSymmetric() const {
	return (symmetric);
}

template<class T>
inline void Matrix<T>::setSymmetric(bool symmetric) {
	this->symmetric = symmetric;
}

template<class T>
Matrix<T>::~Matrix()
{
  for(unsigned int i = 0; i < nRow; i++)
    delete memory[i];
  delete[] memory;
}

template<class T>
void Matrix<T>::reset() {
	memset(memory, 0.0, (nCol * nRow) * sizeof(T));
}

template<class T>
void Matrix<T>::copy(Matrix<T>& a) {
	this->nCol = a.nC();
	this->nRow = a.nR();
	memory = new T*[a.nR()];
	for(unsigned int i = 0; i < a.nR(); i++)
	  memory[i] = new T[a.nC()];

	memcpy(memory,a.memory,sizeof(T)*a.nC()*a.nR());
	symmetric = a.symmetric;
}

template<class T>
ostream& operator<<(ostream &out, Matrix<T> &M) {
	// Since operator<< is a friend of the Point class, we can access
	// Point's members directly.
	for (unsigned int i = 0; i < M.nR(); ++i) {
		for (unsigned int j = 0; j < M.nC(); j++) {
			out << M(i, j) << M.del;
		}
		out << endl;
	}
	return (out);
}


#endif /* MATRIX_H_ */