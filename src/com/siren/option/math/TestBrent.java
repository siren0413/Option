package com.siren.option.math;

import java.util.Arrays;
import java.util.Random;

public class TestBrent {

	public static void main(String[] args) throws Exception {
		System.out.println("haa:" + Black76.price("c", 30, 32, 1, 0.05, 0.3));
		System.out.println("vol:" + Black76.volApproxilation(1.0, 30, 0.05, 1));
		System.out.println("trueVol:" + Black76.impliedVol2(0.3, "c", 30, 32, 1, 0.05));
		System.out.println("n(0.1295)="+Distribution.n(0.1295));
		testBlack76ImpliedVol();

	}

	public static void testBlack76ImpliedVol() throws Exception {

		double s = 30;
		double x = 32;
		double t = 1;
		double r = 0.05;
		double price = 3.22;
		Random random = new Random(System.currentTimeMillis());
		System.out.println("***********");
		// System.out.println(Black76.impliedVol(price, "c", s, x, t, r));
		 System.out.println(Black76.impliedVol2(price, "c", s, x, t, r));
		// System.out.println(Black76.impliedVol3(price, "c", s, x, t, r));
		// System.out.println(Black76.impliedVol4(price, "c", s, x, t, r));
		// System.out.println(Black76.impliedVol5(price, "c", s, x, t, r));
		 System.out.println(Black76.impliedVol6(price, "c", s, x, t, r));
		System.out.println("***********");

		double[] ss = new double[1000000];
		for (int i = 0; i < ss.length; i++) {
			ss[i] = Math.abs(s + (random.nextDouble() - 0.5) * 8.0);
		}

		// System.out.println(Arrays.toString(ss));

		long start;

		// bisection
		System.out.println("Bisection:");
		start = System.currentTimeMillis();
		for (int i = 0; i < ss.length; i++) {
			double result = Black76.impliedVol2(price, "c", ss[i], x, t, r);
			//System.out.println(result);
		}
		System.out.println("Bisection finish:" + (System.currentTimeMillis() - start));

		// Ridder
		System.out.println("Ridder:");
		start = System.currentTimeMillis();
		for (int i = 0; i < ss.length; i++) {
			double result = Black76.impliedVol4(price, "c", ss[i], x, t, r);
			// System.out.println(result);
		}
		System.out.println("Ridder finish:" + (System.currentTimeMillis() - start));

		// Newton
		System.out.println("Newton:");
		start = System.currentTimeMillis();
		for (int i = 0; i < ss.length; i++) {
			double result = Black76.impliedVol6(price, "c", ss[i], x, t, r);
//			if( Math.abs(result - Black76.impliedVol2(price, "c", ss[i], x, t, r)) >0.1) {
//				System.err.println("Newton:"+result+" bisection:"+Black76.impliedVol2(price, "c", ss[i], x, t, r));
//			}
			// System.out.println(result);
		}
		System.out.println("Newton finish:" + (System.currentTimeMillis() - start));
		
		// Combile Newton && bisection
		System.out.println("Newton-Combile:");
		start = System.currentTimeMillis();
		for (int i = 0; i < ss.length; i++) {
			double result = Black76.impliedVol7(price, "c", ss[i], x, t, r);
			// System.out.println(result);
		}
		System.out.println("Newton-Combile finish:" + (System.currentTimeMillis() - start));

		// Secant
		System.out.println("Secant:");
		start = System.currentTimeMillis();
		for (int i = 0; i < ss.length; i++) {
			// double result = Black76.impliedVol5(price, "c", ss[i], x, t, r);
			// System.out.println(result);
		}
		System.out.println("Secant finish:" + (System.currentTimeMillis() - start));

		// brent
		System.out.println("Brent:");
		start = System.currentTimeMillis();
		for (int i = 0; i < ss.length; i++) {
			Black76.impliedVol(price, "c", ss[i], x, t, r);
		}
		System.out.println("Brent finish:" + (System.currentTimeMillis() - start));

		// combile
		System.out.println("Combile:");
		start = System.currentTimeMillis();
		for (int i = 0; i < ss.length; i++) {
			Black76.impliedVol3(price, "c", ss[i], x, t, r);
		}
		System.out.println("Combile finish:" + (System.currentTimeMillis() - start));

	}

	public static double testSin(double x) {
		return Math.sin(x);
	}

}
