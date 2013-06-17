package com.siren.option.math;

public class Distribution {
	
	
	
	static final double p0 = 220.2068679123761;
	static final double p1 = 221.2135961699311;
	static final double p2 = 112.0792914978709;
	static final double p3 = 33.91286607838300;
	static final double p4 = 6.373962203531650;
	static final double p5 = .7003830644436881;
	static final double p6 = .3526249659989109E-01;

	static final double q0 = 440.4137358247522;
	static final double q1 = 793.8265125199484;
	static final double q2 = 637.3336333788311;
	static final double q3 = 296.5642487796737;
	static final double q4 = 86.78073220294608;
	static final double q5 = 16.06417757920695;
	static final double q6 = 1.755667163182642;
	static final double q7 = .8838834764831844E-1;

	static final double cutoff = 7.071;
	static final double root2pi = 2.506628274631001;
	

	public static double CND(double z) {
		double zabs = 0.0;
		double p;
		double expntl, pdf;


		zabs = Math.abs(z);


		// |z| > 37
		if (z > 37.0) {

			p = 1.0;

			return p;

		}

		if (z < -37.0) {

			p = 0.0;

			return p;

		}

		// |z| <= 37.

		expntl = Math.exp(-.5 * zabs * zabs);

		pdf = expntl / root2pi;

		// |z| < cutoff = 10/sqrt(2).

		if (zabs < cutoff) {

			p = expntl * ((((((p6 * zabs + p5) * zabs + p4) * zabs + p3) * zabs + p2) * zabs + p1) * zabs + p0)
					/ (((((((q7 * zabs + q6) * zabs + q5) * zabs + q4) * zabs + q3) * zabs + q2) * zabs + q1) * zabs + q0);

		} else {

			p = pdf / (zabs + 1.0 / (zabs + 2.0 / (zabs + 3.0 / (zabs + 4.0 / (zabs + 0.65)))));

		}

		if (z < 0.0) {

			return p;

		} else {

			p = 1.0 - p;

			return p;

		}

	}
	
	public static double n(double d1) {
		return (1/Math.sqrt(2*Math.PI)) * Math.exp(-d1*d1/2);
	}


}
