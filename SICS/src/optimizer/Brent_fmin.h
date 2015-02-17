#ifndef OPTIMIZER_BRENT_FMIN_
#define OPTIMIZER_BRENT_FMIN_

#include <math.h>
#include <float.h>
#include <type/Constant.h>
#include <model/Model.h>
#define SIGR (3. - sqrt(5.)) * .5

enum variables{v, w, x, u, fv, fw, fx, fu};
const double * default_interval = new double[2] { -5, 5 };

double static Brent_fmin(double * interval, double tol,
		double (*function)(double, vector<char>, int, Model*),
		vector<char> pattern, int node, Model *model, int maxIters) {

	double d, e, p, q, r, u, v, w, x;
	double t2, fu, fv, fw, fx, xm, tol1, tol3;

	tol1 = DBL_EPSILON; + 1.;/* the smallest 1.000... > 1 */

	v = interval[0] + SIGR * (interval[1] - interval[0]);
	w = v;
	x = v;

	d = 0.;/* -Wall */
	e = 0.;
	fx = (*function)(x, pattern, node, model);
	fv = fx;
	fw = fx;
	tol3 = tol / 3.;

	/*  main loop starts here ----------------------------------- */
	int iter = 0;
	while (true) {
		iter++;
		xm = (interval[0] + interval[1]) * .5;
		tol1 = sqrt(DBL_EPSILON) * fabs(x) + tol3;
		t2 = tol1 * 2.;

		/* check stopping criterion */

		if (fabs(x - xm) <= t2 - (interval[1] - interval[0]) * .5) // || iter > maxIters)
			break;
		if (fabs(e) > tol1) { /* fit parabola */

			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = (q - r) * 2.;
			if (q > 0.)
				p = -p;
			else
				q = -q;
			r = e;
			e = d;
		}

		if (fabs(p) >= fabs(q * .5 * r) || p <= q * (interval[0] - x)
				|| p >= q * (interval[1] - x)) { /* a golden-section step */

			if (x < xm)
				e = interval[1] - x;
			else
				e = interval[0] - x;
			d = SIGR * e;
		} else { /* a parabolic-interpolation step */

			d = p / q;
			u = x + d;

			/* f must not be evaluated too close to ax or bx */

			if (u - interval[0] < t2 || interval[1] - u < t2) {
				d = tol1;
				if (x >= xm)
					d = -d;
			}
		}

		/* f must not be evaluated too close to x */

		if (fabs(d) >= tol1)
			u = x + d;
		else if (d > 0.)
			u = x + tol1;
		else
			u = x - tol1;

		fu = (*function)(u, pattern, node, model);

		/*  update  a, b, v, w, and x */

		if (fu <= fx) {
			if (u < x)
				interval[1] = x;
			else
				interval[0] = x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		} else {
			if (u < x)
				interval[0] = u;
			else
				interval[1] = u;
			if (fu <= fw || w == x) {
				v = w;
				fv = fw;
				w = u;
				fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}

	/* end of main loop */
	return (x);
}

#endif /* OPTIMIZER_BRENT_FMIN_ */
