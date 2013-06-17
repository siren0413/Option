package com.siren.option.math;

public class Black76 {

	private static final int ITMAX = 100;
	private static final double EPSILON = 0.00001;

	public static double price(String cp, double s, double x, double t, double r, double v) {

		final double d1 = ((Math.log(s / x)) + (Math.pow(v, 2) / 2) * t) / (v * Math.sqrt(t));
		final double d2 = d1 - v * Math.sqrt(t);
		if ("c".equals(cp)) {
			return Math.exp(-r * t) * (s * Distribution.CND(d1) - x * Distribution.CND(d2));
		} else if ("p".equals(cp)) {
			return Math.exp(-r * t) * (x * Distribution.CND(-d2) - s * Distribution.CND(-d1));
		}
		return -1.0;
	}

	// brent
	public static double impliedVol(double p0, String cp, double S, double X, double T, double R) throws Exception {
		final double price = p0;
		final double tol = 0.0001;
		double a = 0.0001;
		double b = 5.0;
		double c = 5.0;
		double d = 0.0;
		double e = 0.0;
		double fa = price(cp, S, X, T, R, a);
		double fb = price(cp, S, X, T, R, b);
		double fc, p, q, r, s, tol1, xm;
		if ((fa > price && fb > price) || (fa < price && fb < price)) {
			throw new RuntimeException("Root must be bracketed in zbrent: fa=" + fa + " fb=" + fb);
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
			fb = price(cp, S, X, T, R, b);
			// System.err.println("brent:"+i);
		}
		return a;
	}

	public static double sign(double a, double b) {
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
	}

	// bisection
	public static double impliedVol2(double p0, String cp, double S, double X, double T, double R) throws Exception {
		double vol = 0.0;
		double low = 0.0000001;
		double high = 5.0;
		double price = Double.MAX_VALUE;
		int count = 0;
		while (Math.abs(p0 - price) > EPSILON) {
			vol = (high + low) / 2;
			price = Black76.price("c", S, X, T, R, vol);
			if (price > p0) {
				high = vol;
			} else {
				low = vol;
			}
			count++;
			if (count == ITMAX) {
				System.err.println("bisection maxout");
				break;
			}
		}
		// System.err.println("bisection:"+count);
		return vol;
	}

	// combile bisection/brent
	public static double impliedVol3(double p0, String cp, double S, double X, double T, double R) throws Exception {
		double vol = 0.0;
		double low = 0.0001;
		double high = 5.0;
		double price;
		int count = 0;
		while (high - low > 1.0) {
			vol = (high + low) / 2;
			price = Black76.price("c", S, X, T, R, vol);
			if (price > p0) {
				high = vol;
			} else {
				low = vol;
			}
			count++;
		}

		price = p0;
		final double tol = 0.0001;
		double a = low;
		double b = high;
		double c = high;
		double d = 0.0;
		double e = 0.0;
		double fa = price(cp, S, X, T, R, a);
		double fb = price(cp, S, X, T, R, b);
		double fc, p, q, r, s, tol1, xm;
		if ((fa > price && fb > price) || (fa < price && fb < price)) {
			throw new RuntimeException("Root must be bracketed in zbrent: fa=" + fa + " fb=" + fb);
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
			fb = price(cp, S, X, T, R, b);
			// System.err.println("brent:"+i);
		}
		return a;
	}

	// Ridder's method
	public static double impliedVol4(double p0, String cp, double S, double X, double T, double R) throws Exception {
		final double x1 = 0.00001;
		final double x2 = 5.0;
		final double xacc = 0.00001;
		double f1 = price("c", S, X, T, R, x1) - p0;
		double f2 = price("c", S, X, T, R, x2) - p0;
		// System.out.println("f1:"+f1+" f2:"+f2);
		if ((f1 > 0.0 && f2 < 0.0) || (f1 < 0.0 && f2 > 0.0)) {
			double xl = x1;
			double xh = x2;
			double ans = Double.MIN_VALUE;
			for (int i = 0; i < ITMAX; i++) {
				final double xm = 0.5 * (xl + xh);
				final double fm = price("c", S, X, T, R, xm) - p0;
				final double s = Math.sqrt(fm * fm - f1 * f2);
				if (s == 0.0)
					return ans;
				final double xnew = xm + (xm - xl) * ((f1 >= f2 ? 1.0 : -1.0) * fm / s);
				if (Math.abs(xnew - ans) <= xacc)
					return ans;
				ans = xnew;
				final double fnew = price("c", S, X, T, R, ans) - p0;
				if (fnew == 0.0)
					return ans;
				if (sign(fm, fnew) != fm) {
					xl = xm;
					f1 = fm;
					xh = ans;
					f2 = fnew;
				} else if (sign(f1, fnew) != f1) {
					xh = ans;
					f2 = fnew;
				} else if (sign(f2, fnew) != f2) {
					xl = ans;
					f1 = fnew;
				} else {
					System.err.println("never get here");
				}
				if (Math.abs(xh - xl) <= xacc)
					return ans;
			}
			System.err.println("exceed max iterations");
		} else {
			if (f1 == 0.0)
				return x1;
			if (f2 == 0.0)
				return x2;
			System.err.println("root must be bracketed x1:" + x1);
		}
		return -1;
	}

	// secant's method
	public static double impliedVol5(double p0, String cp, double S, double X, double T, double R) throws Exception {
		double xl, rts;
		double x1 = 0.00001;
		double x2 = 5.0;
		double xacc = 0.00001;
		double fl = price("c", S, X, T, R, x1) - p0;
		double f = price("c", S, X, T, R, x2) - p0;
		if (Math.abs(fl) < Math.abs(f)) {
			rts = x1;
			xl = x2;
			swap(fl, f);
		} else {
			xl = x1;
			rts = x2;
		}
		for (int i = 0; i < ITMAX; i++) {
			double dx = (xl - rts) * f / (f - fl);
			xl = rts;
			fl = f;
			rts += dx;
			f = price("c", S, X, T, R, rts) - p0;
			if (Math.abs(dx) < xacc || f == 0.0)
				return rts;
		}
		System.err.println("secant exceed max iteratons");

		return -1;
		/**
		 * // Local variables double x0 = 0.0000001; double x1 = 5.0; double tol
		 * = 0.01; double x, // Calculated value of x at each iteration f0, //
		 * Function value at x0 f1, // Function value at x1 fx, // Function
		 * value at calculated value of x root; // Root, if within desired
		 * tolerance // Set initial function values f0 = price("c", S, X, T, R,
		 * x0) - p0; f1 = price("c", S, X, T, R, x1) - p0; // Loop for finding
		 * root using Secant Method for (int i = 0; i < ITMAX; i++) { x = x1 -
		 * f1 * ((x1 - x0) / (f1 - f0)); fx = price("c", S, X, T, R, x) - p0; x0
		 * = x1; x1 = x; f0 = f1; f1 = fx; // Check whether calculated value is
		 * within tolerance if (Math.abs(x1 - x0) < tol) { root = x1; return
		 * root; } // end if } // end for return x1;
		 */

	}

	// newton's method
	public static double impliedVol6(double p0, String cp, double S, double X, double T, double R) throws Exception {
		double vi = Math.sqrt(Math.abs(Math.log(S / X)) * 2 / T);
		double ci = price("c", S, X, T, R, vi);
		double vega = vega(S, X, T, R, vi);
		double minDiff = Math.abs(p0 - ci);
		int count = 0;
		while (Math.abs(p0 - ci) >= EPSILON && Math.abs(p0 - ci) <= minDiff) {
			count++;
			vi = vi - (ci - p0) / vega;
			ci = price("c", S, X, T, R, vi);
			vega = vega(S, X, T, R, vi);
			minDiff = Math.abs(p0 - ci);
		}
		// System.out.println(count);
		return vi;

	}

	// combile newton && bisection
	public static double impliedVol7(double p0, String cp, double S, double X, double T, double R) throws Exception {

		final double p1 = Black76.price("c", S, X, T, R, 0.1);
		final double p2 = Black76.price("c", S, X, T, R, 0.8);
		if (p1 < p0 && p2 > p0) {
			double vi = Math.sqrt(Math.abs(Math.log(S / X)) * 2 / T);
			double ci = price("c", S, X, T, R, vi);
			double vega = vega(S, X, T, R, vi);
			double minDiff = Math.abs(p0 - ci);
			while (Math.abs(p0 - ci) >= EPSILON && Math.abs(p0 - ci) <= minDiff) {
				vi = vi - (ci - p0) / vega;
				ci = price("c", S, X, T, R, vi);
				vega = vega(S, X, T, R, vi);
				minDiff = Math.abs(p0 - ci);
			}
			// System.out.println(count);
			return vi;
		}

		double vol = 0.0;
		double low = 0.0000001;
		double high = 5.0;
		double price = Double.MAX_VALUE;
		int count = 0;
		while (Math.abs(p0 - price) > EPSILON) {
			vol = (high + low) / 2;
			price = Black76.price("c", S, X, T, R, vol);
			if (price > p0) {
				high = vol;
			} else {
				low = vol;
			}
			count++;
			if (count == ITMAX) {
				System.err.println("bisection maxout");
				break;
			}
		}
		// System.err.println("bisection:"+count);
		return vol;

	}

	public static void swap(double a, double b) {
		double temp = a;
		a = b;
		b = temp;
	}

	public static double volApproxilation(double price, double s, double r, double t) {
		return (price * Math.sqrt(2 * Math.PI)) / (s * Math.exp(-r * t) * Math.sqrt(t));
	}

	public static double vega(double s, double x, double t, double r, double v) {
		final double delay1 = price("c", s, x, t, r, v);
		final double delay2 = price("c", s, x, t, r, v + 0.01);
		final double d1 = ((Math.log(s / x)) + (Math.pow(v, 2) / 2) * t) / (v * Math.sqrt(t));
		return s * Math.exp(-r * t) * Distribution.n(d1) * Math.sqrt(t);
	}

}
