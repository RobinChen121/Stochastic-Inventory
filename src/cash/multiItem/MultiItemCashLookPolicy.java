/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 18, 2019, 9:49:16 PM
 * @Desc: This class is to build a stochastic dynamic programming model for a cash constrained problem 
 *        with two products
 *        
 *        looking for the optimal policy by varying initial cash 
 *
 *
 * 
 */
package cash.multiItem;

import java.util.ArrayList;
import java.util.function.Function;

import sdp.cash.multiItem.Actions;
import sdp.cash.multiItem.CashRecursionMulti;
import sdp.cash.multiItem.CashSimulationMulti;
import sdp.cash.multiItem.CashStateMulti;
import sdp.cash.multiItem.Demands;
import sdp.cash.multiItem.GetPmfMulti;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToExcel;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.probdistmulti.BiNormalDist;


public class MultiItemCashLookPolicy {


	public static void main(String[] args) {
		double[] price = {10, 5};
		double[] variCost = {4, 2};  // higher margin vs lower margin
		
		double iniCash = 25;  // initial cash
		int iniInventory1 = 0;  // initial inventory
		int iniInventory2 = 0;
		
		
		double[][] demand = {{ 6, 6}, {8, 8}}; // higher average demand vs lower average demand
		double[] coe = {0.5, 0.25}; // higher variance vs lower variance
		
		
		double[] salPrice = {2, 1};
		
		int T = demand[0].length; // horizon length
		
		double truncationQuantile = 0.999;
		int stepSize = 1;
		double minCashState = 0;
		double maxCashState = 10000;
		int minInventoryState = 0;
		
		
		int maxInventoryState = 200;
		int Qbound = 100;
		double discountFactor = 1;
		
		
		double Rmin = 25; double Rmax = 80; int incre = 2;
		
		int rowNum = (int) ((Rmax - Rmin)/incre) + 2;
		int row = 0;
		double[][] optResults = new double[rowNum][5];
		
		for (iniCash = Rmin; iniCash <= Rmax; iniCash = iniCash + incre) {
						
		// get demand possibilities for each period
		Distribution[][] distributions = new Distribution[demand.length][T];
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < demand.length; i++){
				distributions[t][i] = new PoissonDist(demand[i][t]);
				
			}
		}
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmfMulti();
		
		// build action list for two items
		Function<CashStateMulti, ArrayList<Actions>> buildActionList = s -> {
			ArrayList<Actions> actions = new ArrayList<>();
			for (int i = 0; i < Qbound; i++)
				for (int j = 0; j < Qbound; j++) {
					if (variCost[0] * i + variCost[1] * j < s.getIniCash() + 0.1) {
						Actions thisAction = new Actions(i, j);
						actions.add(thisAction);
					}					
				}
			return actions;
		};
		
		// Immediate Value Function	      
		ImmediateValueFunction<CashStateMulti, Actions, Demands, Double> immediateValue
		= (IniState, Actions, RandomDemands) -> {
			double action1 = Actions.getFirstAction();
			double action2 = Actions.getSecondAction();
			double demand1 = RandomDemands.getFirstDemand();
			double demand2 = RandomDemands.getSecondDemand();
			double endInventory1 = Math.max(0, IniState.getIniInventory1() + action1 - demand1);
			double endInventory2 = Math.max(0, IniState.getIniInventory2() + action2 - demand2);
			double revenue = price[0] * (IniState.getIniInventory1() + action1 - endInventory1)
					+ price[1] * (IniState.getIniInventory2() + action2 - endInventory2);
			double orderingCosts = variCost[0] * action1 + variCost[1] * action2;
			double salValue = 0;
			if (IniState.getPeriod() == T - 1) {
				salValue = salPrice[0] * endInventory1 + salPrice[1] * endInventory2;
			}
			return revenue - orderingCosts + salValue;
		};
	    	
		// State Transition Function

		StateTransitionFunction<CashStateMulti, Actions, Demands, CashStateMulti> stateTransition = (IniState, Actions, RandomDemands) -> {
			int endInventory1 = IniState.getIniInventory1() + Actions.getFirstAction() - RandomDemands.getFirstDemand();
			endInventory1 = Math.max(0, endInventory1);
			int endInventory2 = IniState.getIniInventory2() + Actions.getSecondAction() - RandomDemands.getSecondDemand();
			endInventory2 = Math.max(0, endInventory2);
			double nextCash = IniState.getIniCash() + immediateValue.apply(IniState, Actions, RandomDemands);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			endInventory1 = endInventory1 > maxInventoryState ? maxInventoryState : endInventory1;
			endInventory2 = endInventory2 < minInventoryState ? minInventoryState : endInventory2;
			nextCash = (int) nextCash;  // rounding states to save computing time
			endInventory1 = (int) endInventory1;
			endInventory2 = (int) endInventory2;
			return new CashStateMulti(IniState.getPeriod() + 1, endInventory1, endInventory2, nextCash);
		};
		
		
		
		/*******************************************************************
		 * Solve
		 */
		CashRecursionMulti recursion = new CashRecursionMulti(discountFactor, pmf, buildActionList,
				                             stateTransition, immediateValue, T);
		int period = 1;
		CashStateMulti iniState = new CashStateMulti(period, iniInventory1, iniInventory2, iniCash);
		long currTime = System.currentTimeMillis();
		double finalValue = iniCash + recursion.getExpectedValueMulti(iniState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first priod is :  Q1 = " + recursion.getAction(iniState).getFirstAction()
				                      + ", Q2 = " + recursion.getAction(iniState).getSecondAction());
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		
		
		/*******************************************************************
		 * Simulating sdp results
		 * 
		 * simulating results a little lower than SDP
		 */
		int sampleNum = 10000;		
		CashSimulationMulti simuation = new CashSimulationMulti(sampleNum, distributions, discountFactor, 
				 recursion, stateTransition, immediateValue);
		double simFinalValue = simuation.simulateSDPGivenSamplNumMulti(iniState);
		System.out.println(simFinalValue);
		
		
		/*******************************************************************
		 * try to find some ordering patters from optTable
		 * 
		 * output results to excel
		 */
		double Q1 = recursion.getAction(iniState).getFirstAction();
		double Q2 = recursion.getAction(iniState).getSecondAction();
		System.out.println("");
		optResults[row][0] = iniInventory1; optResults[row][1] = iniInventory2;
		optResults[row][2] = iniCash; optResults[row][3] = Q1; optResults[row][4] = Q2;
		row++;
		}
		System.out.println("**************************************************");
		
		WriteToExcel wr = new WriteToExcel();
		String fileName = "optTable2.xls";
		String headString =  "x1" + "\t" + "x2" + "\t" + "R" + "\t" + "Q1"+ "\t" + "Q2";
		wr.writeArrayToExcel(optResults, fileName, headString);
		
	}

}
