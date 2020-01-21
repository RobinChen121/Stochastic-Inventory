package capacitated.fitss;

import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.IntStream;

import milp.MIPFitsS;
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
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 *@author: Zhen Chen
 *@email: 15011074486@163.com
 *@date: Jul 14, 2018---1:06:15 PM
 *@description:  fit (s, S) levels for capacitated lot sizing problem
 *
 * may be not k-convex
 * 
 * computation time for 20 periods may be 4s
 * 
 */

public class LevelFitsS {
	public static void main(String[] args) {
		double truncationQuantile = 0.9999;
		double stepSize = 1;
		double minInventory = -150;
		double maxInventory = 400;

		double fixedOrderingCost = 100;
		double variOrderingCost = 0;
		double penaltyCost = 10;
		double[] meanDemand = {6, 6, 6, 6};
		double holdingCost = 1;
		int maxOrderQuantity = 200;
		double coeValue = 0.5;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				//.mapToObj(i -> new PoissonDist(meanDemand[i]))
                .mapToObj(i -> new NormalDist(meanDemand[i], coeValue * meanDemand[i]))
				.toArray(Distribution[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		
//		int T = 8;
//		double[] values = {6, 7};
//		double[] probs = {0.95, 0.05};
//		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
//		.mapToObj(i -> new DiscreteDistribution(values, probs, values.length)) // can be changed to other distributions
//		.toArray(DiscreteDistribution[]::new);	
//		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		
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
		double opt = recursion.getExpectedValue(initialState);
		System.out.println("final optimal expected value is: " + opt);
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
		MIPFitsS findsS = new MIPFitsS(maxOrderQuantity, T);
		double[][] optsS = findsS.getSinglesS(optTable);
		System.out.println("single s, S level: " + Arrays.deepToString(optsS));
		optsS = findsS.getTwosS(optTable);
		System.out.println("two s, S level: " + Arrays.deepToString(optsS));
		optsS = findsS.getThreesS(optTable);
		System.out.println("three s, S level: " + Arrays.deepToString(optsS));
		double sim1 = simuation.simulateSinglesS(initialState, optsS, maxOrderQuantity);
		System.out.printf("one level gap is %.2f%%\n", (sim1 - opt)*100/opt);
		double sim2 = simuation.simulateTwosS(initialState, optsS, maxOrderQuantity);
		System.out.printf("two level gap is %.2f%%\n", (sim2 - opt)*100/opt);
		double sim3 = simuation.simulateThreesS(initialState, optsS, maxOrderQuantity);
		System.out.printf("three level gap is %.2f%%\n", (sim3 - opt)*100/opt);
	}
}

