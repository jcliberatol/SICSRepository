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
	int nCol;
	int nRow;
	T m(char);
	T get3x3determinant ();

public:
	//bool transposed;
	bool symmetric;
	T *memory;
	int ld;
	static char del;
	Matrix(T**, int, int);/*Importing a 2D array*/
	Matrix(); /**Empty object*/
	Matrix(int, int); /**Two dimensional Matrix Constructor allocates memory*/
	Matrix(Matrix<T>&); /**Copy constructor*/
	Matrix(char I, int size); /**Create special kinds of matrices (dense identity)*/
	void reset();/**Reset method, puts all entries in zeros*/
	void transpose ();/**Transposes the matrix, notice it does not perform memory transpose, only index transpose*/
	void copy(Matrix<T>&);/**Copy constructor*/
	int nR(); /** Returns number of rows */
	int nC(); /** Returns number of columns */
	void setIndex(int, int, T);
	T getIndex(int, int);
	T sum(); /** Returns the sum of all objects */
	inline T & operator()(const int nCol, const int nRow); /** Accessing operator for a element */
	T & operator()(const int element); /**Accessing operator for a element */
	friend ostream& operator<<<T>(ostream &, Matrix<T> &); /** Output operator */
	bool isSymmetric() const; /** Symmetry flag for optimizations */
	void setSymmetric(bool symmetric); /** Set to true the symmetry flag */
	virtual ~Matrix();
};

template<class T>
char Matrix<T>::del = ' ';

template<class T>
T Matrix<T>::m(char c){
	return (memory[c - 'a']);
}

template<class T>
T Matrix<T>::sum() {
	T sum = 0;

	for (int i = 0; i < nRow * nCol; i++) {
		sum += memory[i];
	}

	return (sum);
}

template<class T>
int Matrix<T>::nR() {
	return (nRow);
}

template<class T>
int Matrix<T>::nC() {
	return (nCol);
}

template<class T>
inline void Matrix<T>::setIndex(int r, int c, T value) {
		memory[nCol * r + c] = value;
}

template<class T>
inline T Matrix<T>::getIndex(int r, int c) {
		return memory[nCol * r + c];
}

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
Matrix<T>::Matrix(int r, int c) {
	nCol = c;
	nRow = r;
	//transposed = false;
	memory = new T[c * r];
	symmetric = false;
	ld = c;
}

template<class T>
Matrix<T>::Matrix(T** mem , int r, int c) {
	nCol = c;
	nRow = r;
	//transposed = false;
	memory = new T[c * r];
	symmetric = false;
	ld = c;

	for (int var = 0; var < r; ++var) {
		for (int j = 0; j < c; ++j) {
			memory[nCol * var + j] = mem[var][j];
		}
	}
}

template<class T>
Matrix<T>::Matrix(char I, int c) {
	nCol = c;
	nRow = c;
	//transposed = false;
	memory = new T[c * c];
	ld = c;
	symmetric = false;
	this->reset();
	if(I=='I'){
		for(int i = 0 ; i < c ; i++){
			(*this)(i,i)=static_cast<T>(1);
		}
	}
	if(I=='R'){
		srand(time(NULL));
		for(int i = 0 ; i < c*c ; i++){
			memory[i] = (T)(rand())/(T)(INT_MAX);
		}
	}
}

template<class T>
inline T & Matrix<T>::operator()(const int r, const int c) {
		return (memory[nCol * r + c]);

}

template<class T>
T & Matrix<T>::operator()(const int el) {
		return (memory[el]);
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
Matrix<T>::~Matrix() { delete[] memory; }

template<class T>
void Matrix<T>::reset() {
	memset(memory, 0.0, (nCol * nRow) * sizeof(T));
}

template<class T>
void Matrix<T>::copy(Matrix<T>& a) {
	this->nCol = a.nC();
	this->nRow = a.nR();
	memory = new T[a.nC() * a.nR()];

	memcpy(memory,a.memory,sizeof(T)*a.nC()*a.nR());
	symmetric = a.symmetric;
}

template<class T>
ostream& operator<<(ostream &out, Matrix<T> &M) {
	// Since operator<< is a friend of the Point class, we can access
	// Point's members directly.
	for (int i = 0; i < M.nR(); ++i) {
		for (int j = 0; j < M.nC(); j++) {
			out << M(i, j) << M.del;
		}
		out << endl;
	}
	return (out);
}


#endif /* MATRIX_H_ */