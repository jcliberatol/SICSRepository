/*
 * BFGSOptimizer.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef BFGSOPTIMIZER_H_
#define BFGSOPTIMIZER_H_
#include <type/Constant.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
/**
 * Optimizer, takes a function to a pointer that specifies
 * the function to minimize, and the second parameter specifies a function pointer
 * to the gradient of the function, these pointers must be carefully crafted
 * to fit the arguments of the function
 * use a wrapper to call the functions if your function does not follow the form
 * of the input parameters
 * */
int static bfgs(double (*&fntomin)(double*, double*, int, int),
		void (*&gradient)(double*, double*, int, int, double*), double * args,
		double * pars, int nvars, int npars, int maxiter)
{
	/*
	 //call of function
	 double result;
	 result = (*fntomin)(args, pars, n);
	 */
	//call of gradient
	/*
	 double result[n];
	 alloc memory in result
	 (*gradient) (args,pars,n,result);
	 */

	//Tolerances to end the algorithm are in constant.h
	//Allocate all needed arrays
	bool accpoint, enough;
	double *g, *t, *X, *c, *B;
	int n, *l;
	int count, funcount, gradcount;
	double f, gradproj;
	int i, j, ilast, iter = 0;
	double s, steplength;
	double D1, D2;
	double fmin = Constant::INFINITE;
	//carefull of negative max iterations
	if (maxiter <= 0)
		//MAX_ITER_REACHED;
		return (1);
	
	//allocate l
	l = new int[nvars];
	n = 0;

	for (i = 0; i < nvars; i++)
		l[i] = i;

	n = nvars;

	//Assign memory
	g = new double[n];
	t = new double[n];
	X = new double[n];
	c = new double[n];
	//Assign memory to the triangular lower matrix
	B = new double[n * (n + 1) / 2];      
	//evaluate the function at the initial points
	f = (*fntomin)(args, pars, nvars, npars);
	
	if (!(f < Constant::INFINITE))
		//BAD_INITIAL_VALUES;
		return (2);

	//the optimal point for now is f.
	fmin = f;
	funcount = gradcount = 1;
	//run the gradient
	(*gradient)(args, pars, nvars, npars, g);
	iter++;
	ilast = gradcount;

	//Main loop
	do
	{
		if (ilast == gradcount)
		{
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < i; j++)
					B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] = 0.0;
				B[(((i + 1) * (i + 2) / 2) - (i - i) - 1)] = 1.0;
			}
		}

		for (i = 0; i < n; i++)
		{
			X[i] = args[l[i]];
			c[i] = g[l[i]];
		}

		gradproj = 0.0;

		for (i = 0; i < n; i++)
		{
			s = 0.0;

			for (j = 0; j <= i; j++)
				s -= B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] * g[l[j]];
			for (j = i + 1; j < n; j++)
				s -= B[(((j + 1) * (j + 2) / 2) - (j - i) - 1)] * g[l[j]];
			
			t[i] = s;
			gradproj += s * g[l[i]];
		}

		/* search direction is downhill */
		if (gradproj < 0.0)
		{
			steplength = 1.0;
			accpoint = false;
			do
			{
				count = 0;
				for (i = 0; i < n; i++)
				{
					args[l[i]] = X[i] + steplength * t[i];
					//reltest is usually 10.
					if (Constant::reltest + X[i] == Constant::reltest + args[l[i]])
						/* no change */
						count++;
				}

				if (count < n)
				{
					f = (*fntomin)(args, pars, nvars, npars);
					funcount++;
					//Wolfe condition 1
					// f < f(x_k) + c1 * ak * pk * grad(f)
					// Where f(x_k) is fmin
					// c1 is acctol
					// ak is steplength
					// grad(f) is gradproj
					accpoint = (f < Constant::INFINITE)
							&& (f<= fmin+ gradproj * steplength	* Constant::acctol); //acctol is 1x10-4
					
					if (!accpoint)
						// Reduce the steplength by the factor if wolfe condition is not met.
						steplength *= Constant::stepredn;
				}
			} while (!(count == n || accpoint));
			
			//abstol is 1x10-5, reltol is 1e-8
			enough = (f > Constant::abstol) && fabs(f - fmin) > Constant::reltol * (fabs(fmin) + Constant::reltol);
									
			/* stop if value if small or if relative change is low */
			if (!enough)
			{
				count = n;
				fmin = f;
			}

			/* making progress */
			if (count < n)
			{
				fmin = f;
				(*gradient)(args, pars, nvars, npars, g);
				gradcount++;
				iter++;
				D1 = 0.0;

				for (i = 0; i < n; i++)
				{
					t[i] = steplength * t[i];
					c[i] = g[l[i]] - c[i];
					D1 += t[i] * c[i];
				}

				if (D1 > 0)
				{
					D2 = 0.0;
					for (i = 0; i < n; i++)
					{
						s = 0.0;
						for (j = 0; j <= i; j++)
							s += B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] * c[j];
						for (j = i + 1; j < n; j++)
							s += B[(((j + 1) * (j + 2) / 2) - (j - i) - 1)] * c[j];
						X[i] = s;
						D2 += s * c[i];
					}

					D2 = 1.0 + D2 / D1;

					for (i = 0; i < n; i++)
						for (j = 0; j <= i; j++)
							B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] +=
						    		(D2 * t[i] * t[j] - X[i] * t[j] - t[i] * X[j]) / D1;
				/* D1 < 0 */
				} else
					ilast = gradcount;
			/* no progress */
			} else {
				if (ilast < gradcount)
				{
					count = 0;
					ilast = gradcount;
				}
			}
		/* uphill search */
		} else {
			count = 0;

			if (ilast == gradcount)
				count = n;
			else
				ilast = gradcount;
			/* Resets unless has just been reset */
		}

		if (iter >= maxiter)
			break;
		/* periodic restart */
		if (gradcount - ilast > 2 * n)
			ilast = gradcount;
	} while (count != n || ilast != gradcount);

	delete [] g;
	delete [] t;
	delete [] X;
	delete [] c;
	delete [] B;
	delete [] l;

	if (iter < maxiter)
		//SUCCESS;
		return (0);
	//N_CONVERGENCE;
	return (3);
}

#endif /* BFGSOPTIMIZER_H_ */
