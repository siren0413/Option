package com.siren.option.math;

public class Binomial {

	private final static double EPSILON = 0.00001;
	private final static int ITMAX = 100;
	
	
	public static double price(double s, double x, double t, double r, double v, int steps) {

		final double timeStep = t / steps;
		// double t_sqrt = Math.sqrt(t);
		// double sigma2 = v * v;
		final double d1 = ((Math.log(s / x)) + (Math.pow(v, 2) / 2) * t) / (v * Math.sqrt(t));
		final double d2 = d1 - v * Math.sqrt(t);

		final double R = Math.exp(r * timeStep);
		double B = Math.exp(-r * timeStep);

		final double probUp = peizerPratt2(d2, steps);
		final double probUp1 = peizerPratt2(d1, steps);
		final double probDown = 1 - probUp;
		final double up = R * probUp1 / probUp;
		final double down = (R - probUp * up) / probDown;

		final double[] payOffValues = new double[steps + 1];
		final double[] nodePrices = new double[steps + 1];

		nodePrices[0] = s * Math.pow(up, steps);
		payOffValues[0] = payOff(nodePrices[0], x);

		for (int i = 1; i < steps; i++) {
			nodePrices[i] = nodePrices[i - 1] / up * down;
			payOffValues[i] = payOff(nodePrices[i], x);
		}

		for (int S = steps - 1; S >= 0; S--) {
			for (int i = 0; i <= S; i++) {
				nodePrices[i] = nodePrices[i] / up;
				final double nodepayOff = payOff(nodePrices[i], x);

				payOffValues[i] = (payOffValues[i] * probUp + payOffValues[i + 1] * probDown) * B;
				payOffValues[i] = nodeValue(payOffValues[i], nodepayOff);
			}
		}

		return payOffValues[0];

	}

	private static double peizerPratt2(double z, double n) {
		return 0.5 + Math.signum(z)
				* Math.sqrt(0.25 - 0.25 * Math.exp(-Math.pow(z / (n + 1.0 / 3.0 + 0.1 / (n + 1.0)), 2) * (n + 1.0 / 6.0)));
	}

	private static double nodeValue(double value, double nodePayoff) {
		return Math.max(value, nodePayoff);
	}

	private static double payOff(double s, double x) {
		return Math.max(0.0, s - x);
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
			price = Binomial.price(S, X, T, R, vol, 60);
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

}
