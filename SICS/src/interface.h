#ifndef _INTERFACE_H
#include <iostream>
#include <type/Matrix.h>
#include <boost/dynamic_bitset.hpp>
#include <type/PatternMatrix.h>
#include <model/Model.h>
#include <model/ModelFactory.h>
#include <model/SICSGeneralModel.h>
#include <estimation/classical/EMEstimation.h>
#include <input/Input.h>
#include <time.h>
#include <trace/Trace.h>


#define _INTERFACE_H
	void estimatingParameters(int **, int, int, char *, int , char *, double, int, bool, double *);
	void profilerOut(Trace*, int );
#endif
