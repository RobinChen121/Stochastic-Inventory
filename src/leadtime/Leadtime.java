package leadtime;

import java.util.Map;
import java.util.function.Function;
import java.util.stream.IntStream;

import sdp.inventory.GetPmf;
import sdp.inventory.LeadtimeState;
import sdp.inventory.LeadtimeRecursion;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.probdist.NormalDist;

/**
*@author: zhenchen
*@date: Jul 23, 2023, 10:09:59 AM
*@desp: TODO
*
*/

public class Leadtime {

	public static void main(String[] args) {
		double initialInventory = 0; 
	      double[] meanDemand = {10, 10, 10};
	      
	      double truncationQuantile = 0.9999;  
	      double stepSize = 1; 
	      double minInventory = -150;
	      double maxInventory = 300;
	      int T = meanDemand.length;

	      double fixedOrderingCost = 0; 
	      double variOrderingCost = 1; 
	      double holdingCost = 2;
	      double penaltyCost = 10;
	      
	      int maxOrderQuantity = 100;
	      
	      Distribution[] distributions = IntStream.iterate(0, i -> i + 1)
	                                              .limit(T)
	                                              .mapToObj(i -> new PoissonDist(meanDemand[i]))
	                                              //.mapToObj(i -> new NormalDist(meanDemand[i], 0.25 * meanDemand[i]))
	                                              .toArray(Distribution[]::new); // replace for loop
	      double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
	      
	   // feasible actions
		Function<LeadtimeState, double[]> getFeasibleAction = s -> {
			double[] feasibleActions = new double[(int) (maxOrderQuantity / stepSize) + 1];
			int index = 0;
			for (double i = 0; i <= maxOrderQuantity; i = i + stepSize) {
				feasibleActions[index] = i;
				index++;
			}
			return feasibleActions;
		};

		// LeadtimeState transition function
		StateTransitionFunction<LeadtimeState, Double, Double, LeadtimeState> stateTransition = (LeadtimeState, action, randomDemand) -> {
			double nextPreQ = action;
			double nextInventory = LeadtimeState.getIniInventory() + LeadtimeState.getPreQ() - randomDemand;
			// double nextInventory = LeadtimeState.getIniInventory() + action - randomDemand;
//			nextInventory = nextInventory > maxInventory ? maxInventory : nextInventory;
//			nextInventory = nextInventory < minInventory ? minInventory : nextInventory;
			return new LeadtimeState(LeadtimeState.getPeriod() + 1, nextInventory, action);
		};

		// immediate value
		ImmediateValueFunction<LeadtimeState, Double, Double, Double> immediateValue = (LeadtimeState, action, randomDemand) -> {
			double fixedCost = 0, variableCost = 0, inventoryLevel = 0, holdingCosts = 0, penaltyCosts = 0;
			fixedCost = action > 0 ? fixedOrderingCost : 0;
			variableCost = variOrderingCost * action;
			inventoryLevel = LeadtimeState.getIniInventory() + LeadtimeState.getPreQ() - randomDemand;
			//inventoryLevel = LeadtimeState.getIniInventory() + action - randomDemand;
			holdingCosts = holdingCost * Math.max(inventoryLevel, 0);
			penaltyCosts = penaltyCost * Math.max(-inventoryLevel, 0);
			double totalCosts = fixedCost + variableCost + holdingCosts + penaltyCosts;
			return totalCosts;
		};
		
	      
		/*******************************************************************
		 * Solve
		 */
		LeadtimeRecursion recursion = new LeadtimeRecursion(pmf, getFeasibleAction, stateTransition, immediateValue);
		int period = 1;
		double iniInventory = 0;
		LeadtimeState initialState = new LeadtimeState(period, iniInventory, 0);
		long currTime = System.currentTimeMillis();
		double opt = recursion.getExpectedValue(initialState);
		System.out.println("final optimal expected value is: " + opt);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		
		double[][] optTable = recursion.getOptTable();	
//		Map<LeadtimeState, Double> cacheValues = recursion.getCacheValues();
		System.out.println("");

	}

}


