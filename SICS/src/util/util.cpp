#include <util/util.h>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

void filelog(const char * file, const  char * message){
	ofstream myfile;
	char str[80];
	strcpy (str,"irtppdebug/");
	strcat (str,file);
	myfile.open (str,  ios::out | ios::app);
	myfile << message;
	myfile.close();
}

void filelog(const char * file, int message){
	ofstream myfile;
	char str[80];
	strcpy (str,"irtppdebug/");
	strcat (str,file);
	myfile.open (str,  ios::out | ios::app);
	myfile << message;
	myfile.close();
}

void filelog(const char * file, double message){
	ofstream myfile;
	char str[80];
	strcpy (str,"irtppdebug/");
	strcat (str,file);
	myfile.open (str,  ios::out | ios::app);
	myfile << message;
	myfile.close();
}
