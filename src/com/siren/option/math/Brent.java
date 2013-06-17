package com.siren.option.math;

import java.lang.reflect.Method;

public class Brent {

	private static final int ITMAX = 100;
	private static final double EPSILON = 0.00001;

	public static double zbrent(Method m, final double x1, final double x2, final double tol, Object obj, Object... args)
			throws Exception {
		final double price = 0.0;
		double a = x1;
		double b = x2;
		double c = x2;
		double d = 0.0;
		double e = 0.0;
		args[0] = a;
		double fa = (double) m.invoke(obj, args);
		args[0] = b;
		double fb = (double) m.invoke(obj, args);
		double fc, p, q, r, s, tol1, xm;
		if ((fa > price && fb > price) || (fa < price && fb < price)) {
			throw new RuntimeException("Root must be bracketed in zbrent");
		}
		fc = fb;
		for (int i = 0; i < ITMAX; i++) {
	
			if ((fb > price && fc > price) || (fb < price && fc < price)) {
				c = a;
				fc = fa;
				e = d = b - a;
			}
			if (Math.abs(fc) < Math.abs(fb)) {
				a = b;
				b = c;
				c = a;
				fa = fb;
				fb = fc;
				fc = fa;

			}
			tol1 = 2.0 * EPSILON * Math.abs(b) + 0.5 * tol;
			xm = 0.5 * (c - b);
			if (Math.abs(xm) <= tol1 || fb == price)
				return b;
			if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
				s = fb / fa;
				if (a == c) {
					p = 2.0 * xm * s;
					q = 1.0 - s;
				} else {
					q = fa / fc;
					r = fb / fc;
					p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
					q = (q - 1.0) * (r - 1.0) * (s - 1.0);
				}
				if (p > 0.0)
					q = -q;
				p = Math.abs(p);
				double min1 = 3.0 * xm * q - Math.abs(tol1 * q);
				double min2 = Math.abs(e * q);
				if (2.0 * p < (min1 < min2 ? min2 : min2)) {
					e = d;
					d = p / q;
				} else {
					d = xm;
					e = d;
				}
			} else {
				d = xm;
				e = d;
			}
			a = b;
			fa = fb;
			if (Math.abs(d) > tol1)
				b += d;
			else
				b += sign(tol1, xm);
			args[0] = b;
			fb = (double) m.invoke(obj, args);
		}

		return a;
	}

	public static double sign(double a, double b) {
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
	}
}
