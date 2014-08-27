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

template<typename T>
class Matrix;

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
	bool transposed;
	bool symmetric;
	T *memory;
	static char del;
	Matrix(); //Empty object
	Matrix(int, int); //Two dimensional Matrix Constructor allocates memory
	Matrix(Matrix<T>&); //Copy constructor
	Matrix(char I, int size); //Create special kinds of matrices (dense identity)
	void reset();
	void transpose ();
	T getDeterminant ();
	int nR(); // Returns number of rows
	int nC(); //Returns number of columns
	T sum(); // Returns the sum of all objects
	T & operator()(const int nCol, const int nRow); //Accessing operator for a element
	friend ostream& operator<<<T>(ostream &, Matrix<T> &); //Output operator
	bool isSymmetric() const;
	void setSymmetric(bool symmetric);
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
Matrix<T>::Matrix() {
	nCol = 0;
	nRow = 0;
	memory = NULL;
	transposed = false;
	symmetric = false;
}

template<class T>
Matrix<T>::Matrix(Matrix<T>& a) {
	nCol = a.nCol;
	nRow = a.nRow;
	memory = new T[nCol * nRow];
	memcpy(memory,a.memory,sizeof(T)*nCol*nRow);
	transposed = false;
	symmetric = false;
}

template<class T>
Matrix<T>::Matrix(int r, int c) {
	nCol = c;
	nRow = r;
	transposed = false;
	memory = new T[c * r];
	symmetric = false;
}

template<class T>
Matrix<T>::Matrix(char I, int c) {
	nCol = c;
	nRow = c;
	transposed = false;
	memory = new T[c * c];
	symmetric = false;
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
T & Matrix<T>::operator()(const int r, const int c) {
	if (!transposed) {
		return (memory[nCol * r + c]);
	}
	else {
		return (memory[nRow * c + r]);
	}
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
Matrix<T>::~Matrix() {
	if (memory != NULL) {
		delete[] memory;
	}
}

template<class T>
void Matrix<T>::reset() {
	memset(memory, 0.0, (nCol * nRow) * sizeof(T));
}

template<class T>
void Matrix<T>::transpose() {
	transposed = !transposed;
	swap (nCol, nRow);
}

template<class T>
T Matrix<T>::getDeterminant () {
	T det = 0.0;

	if ( nCol == 3 && nRow == 3 ) {
		det = get3x3determinant();
	}

	return (det);
}

template<class T>
T Matrix<T>::get3x3determinant () {
	T det;

	det = m('a')*m('e')*m('i') + m('d')*m('h')*m('c') + m('g')*m('b')*m('f');
	det -= m('g')*m('e')*m('c') + m('a')*m('h')*m('f') + m('d')*m('b')*m('i');

	return (det);
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
