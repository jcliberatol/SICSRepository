//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Main.h"

using namespace std;

int main() {
	PatternMatrix m;
	boost::dynamic_bitset<> l(5);
	for (int var = 0; var < 5; ++var) {
		cout<<l;
		m.push(l);
		l[var] = 1;
		cout<<endl;
	}
	cout<<endl;
	cout<<m;
	cout<<endl;
	cout<<m(l);
	cout<<endl;

	m.flush();
	Input input;
	input.importCSV((char *) "input.csv", m, 0, 0 );
	Matrix<double> mat(3,3);
	input.importCSV((char * ) "decimals.csv",mat,0,0);
	cout << mat;
	Matrix<double> cuad(41,2);
	input.importCSV((char * ) "Cuads.csv",cuad,1,0);
	cout << cuad;
	return 0;
}
