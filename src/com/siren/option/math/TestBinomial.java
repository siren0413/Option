package com.siren.option.math;

public class TestBinomial {

	public static void main(String[] args) throws Exception {
		System.out.println("binomial:" + Binomial.price(30, 32, 1, 0.05, 0.3, 60));
		System.out.println("Black76:" + Black76.price("c", 30, 32, 1, 0.05, 0.3));
		
		
		long start = 0;
		System.out.println("bimonial:");
		start = System.currentTimeMillis();
		for(int i = 0 ; i < 100000; i ++) {
			Binomial.price(30, 32, 1, 0.05, 0.3, 60);
		}
		System.out.println("binomial finish:"+ (System.currentTimeMillis() - start));
		
		
		
		System.out.println("Black76:");
		start = System.currentTimeMillis();
		for(int i = 0 ; i < 100000; i ++) {
			Black76.price("c", 30, 32, 1, 0.05, 0.3);
			Black76.impliedVol6(3.22, "c", 30, 32, 1, 0.05);
		}
		System.out.println("Black76 finish:"+ (System.currentTimeMillis() - start));
		
		
		
	}
}
 