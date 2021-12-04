/**
 * @date: Nov 21, 2020
 */
package cash.risk;

import java.text.DecimalFormat;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import cash.strongconstraint.FindsCS.FindCCrieria;
import sdp.cash.CashRecursion;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Nov 21, 2020
 * @Desc: stochastic dynamic programming model to maximize the survival probability, 
 *        the paper is Archibald & Betts (2002), management science.
 * 
 *
 */
public class cashSurvival {

	/**
	 * @param args
	 * @date: Nov 21, 2020, 6:01:10 PM 
	 */
	public static void main(String[] args) {
		double[] meanDemand = {5, 5, 5, 5};
		//double[] meanDemand = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
		double iniInventory = 0;
		double iniCash = 15;
		double fixOrderCost = 0;
		double variCost = 1;
		double[] price = {3, 3, 3, 3};
		double depositeRate = 0;
		double salvageValue = 0.5;
		double holdingCost = 0;	
		FindCCrieria criteria = FindCCrieria.XRELATE;	
		
		double overheadCost = 10; // costs like wages or rents which is required to pay in each period
		double overheadRate = 0; // rate from revenue to pay overhead wages
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash
		
		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -1000; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;
		double penaltyCost = 0; // can also be overdraft rate; large penalty cost cause big gaps for simulation results, since may generate zero demand;
		double discountFactor = 1;		
		
		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				// .mapToObj(i -> new NormalDist(meanDemand[i], Math.sqrt(meanDemand[i]))) // can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i])).toArray(Distribution[]::new);
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		

		// feasible actions
		// in fact, no cash constraint in this paper
		Function<CashState, double[]> getFeasibleAction = s -> {
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxOrderQuantity + 1).toArray();
		};
		
		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			int t = state.getPeriod() - 1;
			double revenue = price[t] * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double deposite = (state.getIniCash() - fixedCost - variableCost) * (1 + depositeRate);
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = (1 - overheadRate) * revenue + deposite - holdCosts - overheadCost
					- state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			double endCash = state.getIniCash() + cashIncrement;
			if (endCash < 0) {
				cashIncrement += penaltyCost * endCash; // can change to overdraft interest
			}
			return cashIncrement;
		};

		// state transition function
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition = (state, action,
				randomDemand) -> {
			double nextInventory = Math.max(0, state.getIniInventory() + action - randomDemand);
			double nextCash = state.getIniCash() + immediateValue.apply(state, action, randomDemand);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			// cash is integer or not
			nextCash = Math.round(nextCash * 1) / 1; // the right should be a decimal
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};
		

		/*******************************************************************
		 * Solve
		 */
		CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition,
				immediateValue, discountFactor);
		int period = 1;		
		CashState initialState = new CashState(period, iniInventory, iniCash);
		long currTime = System.currentTimeMillis();
		recursion.setTreeMapCacheAction();
		double finalValue = recursion.getSurvProb(initialState);
		System.out.println("survival probability for this initial state is: " + finalValue);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		double[][] optTable = recursion.getOptTable();
		System.out.println();
		
		/*******************************************************************
		 * Simulate the result
		 */
		
		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue2 = (state, action, randomDemand) -> {
			int t = state.getPeriod() - 1;
			double revenue = price[t] * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double deposite = (state.getIniCash() - fixedCost - variableCost) * (1 + depositeRate);
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = (1 - overheadRate) * revenue + deposite - holdCosts - overheadCost
					- state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			double endCash = state.getIniCash() + cashIncrement;
			return cashIncrement;
		};
		
		int sampleNum = 100000;
		CashSimulation simulation = new CashSimulation(distributions, sampleNum, recursion, discountFactor); // no need to add overheadCost in this class
		double[] result = simulation.simulateSDPGivenSamplNumLostRate(initialState, immediateValue2);
		DecimalFormat df2 = new DecimalFormat("###, ###");
		System.out.println("\nfinal simulated survival probability in " + df2.format(sampleNum) + " samples is: " + result[0]);
		System.out.println("\nfinal simulated lost sale rate " + " is: " + result[1]);
		 
		
	}

}







