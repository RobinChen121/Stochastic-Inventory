package cash.overdraft;

import java.util.Arrays;

/**
*@author: zhenchen
*@date: Oct 17, 2023, 8:42:51 AM
*@desp: TODO
*
*/

public class MultiProduct {
	public static void main(String[] args) {
		double[] price = {5, 10};
		double[] variCost = {1, 2};  // higher margin vs lower margin
		double depositeRate = 0;
		
		
		double iniCash = 0;  // initial cash
		int iniI1 = 0;  // initial inventory
		int iniI2 = 0;
		
		// gamma distribution:mean demand is shape / beta and variance is shape / beta^2
		// beta = 1 / scale
		// shape = demand * beta
		// variance = demand / beta
		// gamma in ssj library: alpha is alpha, and lambda is beta(beta)
		int T = 4; // horizon length
		double[] meanDemands = new double[] {10, 5};		
		double[][] demand = new double[2][T]; // higher average demand vs lower average demand
		double[] beta = {10, 1}; // lower variance vs higher variance
		
		double v1 = variCost[0]; double v2 = variCost[1];
		double p1 = price[0]; double p2 = price[1];
		for (int t = 0; t < T; t++) {
			demand[0][t] = meanDemands[0];
			demand[1][t] = meanDemands[1];
		}
				
		double[] salValueUnit = Arrays.stream(variCost).map(a -> a*0.5).toArray();
		int m = demand.length; // number of products
			
		double truncationQuantile = 0.9999; // may affect poisson results
		int stepSize = 1;
		double minCashState = 0;
		double maxCashState = 10000;
		int minInventoryState = 0;	
		int maxInventoryState = 200;
		int Qbound = 50;
		double discountFactor = 1;
		
		
	}

}


