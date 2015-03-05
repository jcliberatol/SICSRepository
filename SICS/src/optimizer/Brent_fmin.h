#ifndef OPTIMIZER_BRENT_FMIN_
#define OPTIMIZER_BRENT_FMIN_

#include <math.h>
#include <float.h>
#include <type/Constant.h>
#include <model/Model.h>
#define SIGR (3. - sqrt(5.)) * .5

enum variables {
	v, w, x, u, fv, fw, fx, fu
};
const double * default_interval = new double[2] { -5, 5 };

double static Brent_fmin(double * interval, double tol,
		double (*function)(double, vector<char>, int, Model*),
		vector<char> pattern, int node, Model *model, int maxIters) {

	const double c = (3. - sqrt(5.)) * .5;
	/* Local variables */
	double a, b, d, e, p, q, r, u, v, w, x;
	double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
	/* eps is approximately the square root of the relative machine precision. */
	eps = DBL_EPSILON;
	tol1 = eps + 1.;/* the smallest 1.000... > 1 */
	eps = sqrt(eps);
	a = interval[0];
	b = interval[1];
	v = a + c * (b - a);
	w = v;
	x = v;
	d = 0.;/* -Wall */
	e = 0.;
	fx = (*function)(x, pattern, node, model);
	fv = fx;
	fw = fx;
	tol3 = tol / 3.;
	/* main loop starts here ----------------------------------- */
	for (;;) {
		xm = (a + b) * .5;
		tol1 = eps * fabs(x) + tol3;
		t2 = tol1 * 2.;
		/* check stopping criterion */
		if (fabs(x - xm) <= t2 - (b - a) * .5)
			break;
		p = 0.;
		q = 0.;
		r = 0.;
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
		if (fabs(p) >= fabs(q * .5 * r) || p <= q * (a - x)
				|| p >= q * (b - x)) { /* a golden-section step */
			if (x < xm)
				e = b - x;
			else
				e = a - x;
			d = c * e;
		} else { /* a parabolic-interpolation step */
			d = p / q;
			u = x + d;
			/* f must not be evaluated too close to ax or bx */
			if (u - a < t2 || b - u < t2) {
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
		/* update a, b, v, w, and x */
		if (fu <= fx) {
			if (u < x)
				b = x;
			else
				a = x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		} else {
			if (u < x)
				a = u;
			else
				b = u;
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
	return x;
}

double static Brent_fmin(double * interval, double tol,
		double (*function)(double, bool *, int, int, Model*), bool * pattern,
		int size, int node, Model *model, int maxIters) {

	const double c = (3. - sqrt(5.)) * .5;
	/* Local variables */
	double a, b, d, e, p, q, r, u, v, w, x;
	double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
	/* eps is approximately the square root of the relative machine precision. */
	eps = DBL_EPSILON;
	tol1 = eps + 1.;/* the smallest 1.000... > 1 */
	eps = sqrt(eps);
	a = interval[0];
	b = interval[1];
	v = a + c * (b - a);
	w = v;
	x = v;
	d = 0.;/* -Wall */
	e = 0.;
	fx = (*function)(x, pattern, size, node, model);
	fv = fx;
	fw = fx;
	tol3 = tol / 3.;
	/* main loop starts here ----------------------------------- */
	for (;;) {
		xm = (a + b) * .5;
		tol1 = eps * fabs(x) + tol3;
		t2 = tol1 * 2.;
		/* check stopping criterion */
		if (fabs(x - xm) <= t2 - (b - a) * .5)
			break;
		p = 0.;
		q = 0.;
		r = 0.;
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
		if (fabs(p) >= fabs(q * .5 * r) || p <= q * (a - x)
				|| p >= q * (b - x)) { /* a golden-section step */
			if (x < xm)
				e = b - x;
			else
				e = a - x;
			d = c * e;
		} else { /* a parabolic-interpolation step */
			d = p / q;
			u = x + d;
			/* f must not be evaluated too close to ax or bx */
			if (u - a < t2 || b - u < t2) {
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
		fu = (*function)(u, pattern, size, node, model);
		/* update a, b, v, w, and x */
		if (fu <= fx) {
			if (u < x)
				b = x;
			else
				a = x;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		} else {
			if (u < x)
				a = u;
			else
				b = u;
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
	return x;
}

#endif /* OPTIMIZER_BRENT_FMIN_ */
