/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 18, 2019, 9:49:16 PM
 * @Desc: This class is to build a stochastic dynamic programming model for a cash constrained problem 
 *        with two products
 *        
 *        Roberto use the states of two consecutive periods to compute state transition probability, by 
 *        get the random demand given consecutive states and actions
 *
 *
 * 
 */
package multiItem;

import java.util.ArrayList;
import java.util.function.Function;


import sdp.cash.multiItem.Action;
import sdp.cash.multiItem.CashStateMulti;
import sun.print.resources.serviceui;
import umontreal.ssj.probdistmulti.BiNormalDist;
import umontreal.ssj.probdistmulti.ContinuousDistributionMulti;

public class MultiItemCash {


	public static void main(String[] args) {
		double[] price = {10, 5};
		double[] variCost = {4, 2};
		
		double iniCash = 30;
		
		double[][] demand = {{8, 8, 8, 8}, {5, 5, 5, 5, 5}};		
		double coe = 0.25;
		
		int T = demand[0].length; // horizon length
		int N = demand.length; // item number
		
		double minCashState = 0;
		double maxCashState = 10000;
		double minInventoryState = 0;
		double maxInventoryState = 200;
		double Qbound = 50;
		
		// get demand possibilities for each period
		ContinuousDistributionMulti[] distributions =  new ContinuousDistributionMulti[T];
		for (int t = 0; t < T; t++)
			distributions[t] = new BiNormalDist(demand[0][t], coe * demand[0][t], demand[1][t], coe * demand[1][t], 0);
		
		// build action list for two items
		Function<CashStateMulti, ArrayList<int[]>> buildActionList = s -> {
			ArrayList<int[]> actions = new ArrayList<>();
			for (int i = 0; i < (int) Qbound; i++)
				for (int j = 0; j < (int) Qbound; j++) {
					if (variCost[0] * i + variCost[1] * j < s.getIniCash() + 0.1) {
						int[] thisAction = {i,j};
						actions.add(thisAction);
					}					
				}
			return actions;
		};

		
		
		// Immediate Value Function	      
	    ImmediateValueFunction<CashStateMulti, Integer, Integer, CashStateMulti, Double> immediateValueFunction
	    	= (IniState, action1, action2, FinalState) -> {
	    		double revenue = price[0] * (IniState.getIniInventory1() + action1 - FinalState.getIniInventory1());
	    		double orderingCosts = variCost[0] * action1 + variCost[1] * action2;
	    		return revenue - orderingCosts;
	    	};
		

	}

}
