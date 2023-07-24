package leadtime;

import java.lang.Thread.State;
import java.util.function.Function;
import java.util.stream.IntStream;

import sdp.inventory.GetPmf;
import sdp.inventory.LeadtimeState;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
*@author: zhenchen
*@date: Jul 23, 2023, 10:09:59 AM
*@desp: TODO
*
*/

public class Leadtime {

	public static void main(String[] args) {
		double initialInventory = 0; 
	      double[] meanDemand = {10, 10};
	      
	      double truncationQuantile = 0.9999;  
	      double stepSize = 1; 
	      double minInventory = -150;
	      double maxInventory = 300;
	      int T = meanDemand.length;

	      double fixedOrderingCost = 0; 
	      double proportionalOrderingCost = 1; 
	      double holdingCost = 2;
	      double penaltyCost = 10;
	      
	      int maxOrderQuantity = 200;
	      
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

		// state transition function
		StateTransitionFunction<LeadtimeState, Double, Double, LeadtimeState> stateTransition = (state, action, randomDemand) -> {
			double nextPreQ = action;
					state.getIniInventory() + action - randomDemand;
			nextInventory = nextInventory > maxInventory ? maxInventory : nextInventory;
			nextInventory = nextInventory < minInventory ? minInventory : nextInventory;
			return new State(state.getPeriod() + 1, nextInventory);
		};

		// immediate value
		ImmediateValueFunction<State, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double fixedCost = 0, variableCost = 0, inventoryLevel = 0, holdingCosts = 0, penaltyCosts = 0;
			fixedCost = action > 0 ? fixedOrderingCost : 0;
			variableCost = variOrderingCost * action;
			inventoryLevel = state.getIniInventory() + action - randomDemand;
			holdingCosts = holdingCost * Math.max(inventoryLevel, 0);
			penaltyCosts = penaltyCost * Math.max(-inventoryLevel, 0);
			double totalCosts = fixedCost + variableCost + holdingCosts + penaltyCosts;
			return totalCosts;
		};
	      

	}

}


