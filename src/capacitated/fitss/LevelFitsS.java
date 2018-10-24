package capacitated.fitss;

import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.IntStream;

import sdp.inventory.CheckKConvexity;
import sdp.inventory.GetPmf;
import sdp.inventory.Drawing;
import sdp.inventory.Recursion;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.Recursion.OptDirection;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
 *@author: Zhen Chen
 *@email: 15011074486@163.com
 *@date: Jul 14, 2018---1:06:15 PM
 *@description:  fit (s, S) levels for capacitated lot sizing problem
 *
 * may be not k-convex
 */

public class LevelFitsS {
	public static void main(String[] args) {
		double truncationQuantile = 1;
		double stepSize = 1;
		double minInventory = -500;
		double maxInventory = 500;

		double fixedOrderingCost = 22;
		double variOrderingCost = 1;
		double penaltyCost = 10;
		double[] meanDemand = { 4.9, 18.8, 6.4, 27.9};
		double holdingCost = 3;
		int maxOrderQuantity = 10;  
		
		double discountFactor = 1;

		// get demand possibilities for each period
//		int T = meanDemand.length;
//		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
//				.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other distributions
//				.toArray(PoissonDist[]::new);
//		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		
		int T = 8;
		double[] values = {6, 7};
		double[] probs = {0.95, 0.05};
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
		.mapToObj(i -> new DiscreteDistribution(values, probs, values.length)) // can be changed to other distributions
		.toArray(DiscreteDistribution[]::new);	
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		
		// feasible actions
		Function<State, double[]> getFeasibleAction = s -> {
			double[] feasibleActions = new double[(int) (maxOrderQuantity / stepSize) + 1];
			int index = 0;
			for (double i = 0; i <= maxOrderQuantity; i = i + stepSize) {
				feasibleActions[index] = i;
				index++;
			}
			return feasibleActions;
		};

		// state transition function
		StateTransitionFunction<State, Double, Double, State> stateTransition = (state, action, randomDemand) -> {
			double nextInventory = state.getIniInventory() + action - randomDemand;
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

		/*******************************************************************
		 * Solve
		 */
		Recursion recursion = new Recursion(OptDirection.MIN, pmf, getFeasibleAction, stateTransition, immediateValue);
		int period = 1;
		double iniInventory = 0;
		State initialState = new State(period, iniInventory);
		long currTime = System.currentTimeMillis();
		System.out.println("final optimal expected value is: " + recursion.getExpectedValue(initialState));
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");

		/*******************************************************************
		 * Simulating sdp results
		 */
		int sampleNum = 10000;
		SimulateFitsS simuation = new SimulateFitsS(distributions, sampleNum, recursion);
		simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001; 
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);
		
		/*******************************************************************
		 * Fit (s, S) levels
		 */
		System.out.println("");
		double[][] optTable = recursion.getOptTable();
		
		FindsS findsS = new FindsS(maxOrderQuantity, T);
		double[][] optsS = findsS.getSinglesS(optTable);
		System.out.println("single s, S level: " + Arrays.deepToString(optsS));
		optsS = findsS.getTwosS(optTable);
		System.out.println("two s, S level: " + Arrays.deepToString(optsS));
		optsS = findsS.getThreesS(optTable);
		System.out.println("three s, S level: " + Arrays.deepToString(optsS));
		simuation.simulateSinglesS(initialState, optsS, maxOrderQuantity);
		simuation.simulateTwosS(initialState, optsS, maxOrderQuantity);
		simuation.simulateThreesS(initialState, optsS, maxOrderQuantity);
		
 		/*******************************************************************
		 * Check K-convexity
		 */	
 		int minInventorys = 0;
		int maxInventorys = 100; 
		int xLength = maxInventorys - minInventorys + 1;
 		double[][] yG = new double[xLength][2];
		int index = 0;
		for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
			yG[index][0] = initialInventory;
			yG[index][1] = recursion.getExpectedValue(new State(period, initialInventory));
			index++;
		}
 		CheckKConvexity.checkCK(yG, fixedOrderingCost, maxOrderQuantity);
 		Drawing.drawXC(yG);
 		Drawing.drawGAndsS(yG, fixedOrderingCost);
	}
}

