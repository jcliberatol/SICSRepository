//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "SICS.h"
#include <OpenBlas/cblas.h>
#include <OpenBlas/common.h>
#include <boost/dynamic_bitset.hpp>
#include <types/PatternMatrix.h>

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

	return 0;
}
