//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "SICS.h"

using namespace std;

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	input iobj;
	bool t = iobj.checkFile("/home/mirt/git/testo");
	cout<<"t := "<<t<<endl;
	return 0;
}
