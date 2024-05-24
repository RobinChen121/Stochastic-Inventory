/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 18, 2019, 9:49:16 PM
 * @Desc: This class is to build a stochastic dynamic programming model for a cash constrained problem 
 *        with two products, the states used are x, R and the action is y;
 *        
 *        4 periods running time is 80s;
 *        2 periods running time is 0.5s;
 *        
 *
 * 
 */
package cash.multiItem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;


import sdp.cash.multiItem.CashRecursionMultiXR;
import sdp.cash.multiItem.CashSimulationMulti;
import sdp.cash.multiItem.CashSimulationMultiXR;
import sdp.cash.multiItem.CashStateMultiXR;

import sdp.cash.multiItem.GetPmfMulti;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.probdistmulti.BiNormalDist;


public class MultiItemCashXR {


	public static void main(String[] args) {
		double[] price = {5, 10};
		double[] variCost = {1, 2};  // higher margin vs lower margin
		
		double iniCash = 0;  // initial cash
		int iniInventory1 = 0;  // initial inventory
		int iniInventory2 = 0;
		
		double depositeRate = 0;
		
		// gamma distribution:mean demand is shape / beta and variance is shape / beta^2
		// beta = 1 / scale
		// shape = demand * beta
		// variance = demand / beta
		// gamma in ssj: alpha is alpha, and lambda is beta(beta)
		int T = 2; // horizon length
		double[] meanDemands = new double[] {20, 10};	
		double[][] demand = new double[2][T]; // higher average demand vs lower average demand
		double[] beta = {10, 1}; // higher variance vs lower variance
		
		double[] salPrice = Arrays.stream(variCost).map(a -> a * 0.5).toArray();	
		
		int m = meanDemands.length; // number of products
		double d1 = meanDemands[0];
		double d2 = meanDemands[1];
		for (int t = 0; t < T; t++) {
			demand[0][t] = d1;
			demand[1][t] = d2;
		}
		
		double truncationQuantile = 0.99; // may affect poisson results
		int stepSize = 1;
		double minCashState = 0;
		double maxCashState = 10000;
		int minInventoryState = 0;	
		int maxInventoryState = 200;
		int Qbound = 50;
		double discountFactor = 1;
		
		// get demand possibilities for each period
		//Distribution[][] distributions =  new GammaDist[m][T];
		Distribution[][] distributions =  new PoissonDist[m][T];
		//Distribution[][] distributions =  new NormalDist[m][T];
		for (int i = 0; i < m; i++)
			for (int t = 0; t < T; t++) {
				// distributions[i][t] = new GammaDist(demand[i][t]* beta[i], beta[i]);
				distributions[i][t] = new PoissonDist(demand[i][t]);
				//distributions[i][t]= new NormalDist(demand[i][t], 0.1 * demand[i][t]);
			}
		
		// build action list (y1, y2) for two items
		Function<CashStateMultiXR, ArrayList<double[]>> buildActionList = s -> {
			ArrayList<double[]> actions = new ArrayList<>();
			int miny1 = (int) s.getIniInventory1();
			int miny2 = (int) s.getIniInventory2();
			for (int i = miny1; i < miny1 + Qbound; i++)
				for (int j = miny2; j < miny2 + Qbound; j++) {
				//	if (variCost[0] * i + variCost[1] * j < s.getIniR() + 0.1) {
						double[] thisActions = {i, j};
						actions.add(thisActions);
				//	}					
				}
			return actions;
		};
		
		// Immediate Value Function	      
		ImmediateValueFunction<CashStateMultiXR, double[], double[], Double> immediateValue
		= (IniState, Actions, RandomDemands) -> {
			double action1 = Actions[0];
			double action2 = Actions[1];
			double demand1 = RandomDemands[0];
			double demand2 = RandomDemands[1];
			double endInventory1 = Math.max(0, action1 - demand1);
			double endInventory2 = Math.max(0, action2 - demand2);
			double revenue1 = price[0] * (action1 - endInventory1);
			double revenue2 = price[1] * (action2 - endInventory2);
			double revenue = revenue1 + revenue2;
			double initialCash = IniState.getIniR() - variCost[0] * IniState.getIniInventory1() - variCost[1] * IniState.getIniInventory2();
			double orderingCostY1 = variCost[0] * action1;
			double orderingCostY2 = variCost[1] * action2;
			double orderingCostsY = orderingCostY1 + orderingCostY2;
			double salValue = 0;
			if (IniState.getPeriod() == T) {
				salValue = salPrice[0] * endInventory1 + salPrice[1] * endInventory2;
			}
			return revenue + (1 - depositeRate) * (IniState.getIniR() - orderingCostsY) + salValue - initialCash;
		};
	    	
		// State Transition Function

		StateTransitionFunction<CashStateMultiXR, double[], double[], CashStateMultiXR> stateTransition = (IniState, Actions, RandomDemands) -> {
			double endInventory1 = Actions[0] - RandomDemands[0];
			endInventory1 = Math.max(0, endInventory1);
			double endInventory2 = Actions[1] - RandomDemands[1];
			endInventory2 = Math.max(0, endInventory2);
			double initialCash = IniState.getIniR() - variCost[0] * IniState.getIniInventory1() - variCost[1] * IniState.getIniInventory2();
			double nextCash = initialCash + immediateValue.apply(IniState, Actions, RandomDemands);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			endInventory1 = endInventory1 > maxInventoryState ? maxInventoryState : endInventory1;
			endInventory2 = endInventory2 < minInventoryState ? minInventoryState : endInventory2;
			nextCash = (int) nextCash;  // rounding states to save computing time
			endInventory1 = (int) endInventory1;
			endInventory2 = (int) endInventory2;
			double nextR = (int) nextCash + variCost[0] * endInventory1 + variCost[1] * endInventory2;
			return new CashStateMultiXR(IniState.getPeriod() + 1, endInventory1, endInventory2, nextR);
		};
		
		GetPmfMulti PmfMulti = new GetPmfMulti(distributions, truncationQuantile, stepSize);
		
		/*******************************************************************
		 * Solve
		 */
		CashRecursionMultiXR recursion = new CashRecursionMultiXR(discountFactor, PmfMulti, buildActionList,
				                             stateTransition, immediateValue, T);
		int period = 1;
		CashStateMultiXR iniState = new CashStateMultiXR(period, iniInventory1, iniInventory2, iniCash);
		long currTime = System.currentTimeMillis();
		double finalValue = iniCash + recursion.getExpectedValue(iniState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first priod is :  y1 = " + recursion.getAction(iniState)[0]
				                      + ", y2 = " + recursion.getAction(iniState)[1]);
		double time = (System.currentTimeMillis() - currTime) / 1000.0;
		System.out.println("running time is " + time + "s");
		
//		
//		
//		/*******************************************************************
//		 * Simulating sdp results
//		 * 
//		 * simulating results a little lower than SDP
//		 */
//		int sampleNum = 10000;		
//		CashSimulationMultiXR simuation = new CashSimulationMultiXR(sampleNum, distributions, discountFactor, 
//				 recursion, stateTransition, immediateValue);
//		double simFinalValue = simuation.simulateSDPGivenSamplNum(iniState);
//		System.out.println(simFinalValue);
//		
//		
//		/*******************************************************************
//		 * try to find some ordering patters from optTable
//		 * 
//		 * output results to excel
//		 */
//		System.out.println("");
//		double[][] optTable = recursion.getOptTable(variCost);
//		WriteToExcel wr = new WriteToExcel();
//		String fileName = "optTable" + "_c1=" + variCost[0] + "c2=" + variCost[1] + ".xls";
//		String headString =  "period" + "\t" + "x1" + "\t" + "x2" + "\t" + "w"+ "\t" + "R" + "\t" + "is limited cash and both ordering" + "\t" + "alpha"
//				 				+ "\t" + "y1"+ "\t" + "y2" + "\t" + "c1" + "\t" + "c2";
//		wr.writeArrayToExcel(optTable, fileName, headString);
//		
	}

}
