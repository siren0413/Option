package com.siren.option.math;

public class Binomial {

	public static double price(double s, double x, double t, double r, double v, int steps) {

		double timeStep = t / steps;
		double t_sqrt = Math.sqrt(t);
		double sigma2 = v * v;
		double d1 = ((Math.log(s / x)) + (Math.pow(v, 2) / 2) * t) / (v * Math.sqrt(t));
		double d2 = d1 - v * Math.sqrt(t);

		double R = Math.exp(r * timeStep);
		double B = Math.exp(-r * timeStep);

		double probUp = peizerPratt2(d2, steps);
		double probUp1 = peizerPratt2(d1, steps);
		double probDown = 1 - probUp;
		double up = R * probUp1 / probUp;
		double down = (R - probUp * up) / probDown;

		double[] payOffValues = new double[steps + 1];
		double[] nodePrices = new double[steps + 1];

		nodePrices[0] = s * Math.pow(up, steps);
		payOffValues[0] = payOff(nodePrices[0], x);

		for (int i = 1; i < steps; i++) {
			nodePrices[i] = nodePrices[i - 1] / up * down;
			payOffValues[i] = payOff(nodePrices[i], x);
		}

		for (int S = steps - 1; S >= 0; S--) {
			for (int i = 0; i <= S; i++) {
				nodePrices[i] = nodePrices[i] / up;
				double nodepayOff = payOff(nodePrices[i], x);

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

}
