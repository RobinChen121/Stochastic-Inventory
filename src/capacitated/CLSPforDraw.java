package capacitated;

import java.util.function.Function;
import java.util.stream.IntStream;

import sdp.inventory.CheckKConvexity;
import sdp.inventory.Drawing;
import sdp.inventory.GetPmf;
import sdp.inventory.Recursion;
import sdp.inventory.Simulation;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.Recursion.OptDirection;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 9, 2018---12:48:56 PM
 * @description: solving capacitated lot sizing problem, and drawing some
 *               pictures for analysis; a best example for draw k convex of
 *               poisson demands : double fixedOrderingCost = 500; double
 *               variOrderingCost = 0; double penaltyCost = 10; double[]
 *               meanDemand = {9, 23, 53, 29}; double holdingCost = 2; int
 *               maxOrderQuantity = 100;
 *
 */

public class CLSPforDraw {
	public static void main(String[] args) {
		double truncationQuantile = 0.9999;
		double stepSize = 1;
		double minInventory = -500;
		double maxInventory = 500;

		double fixedOrderingCost = 100;
		double variOrderingCost = 0;
		double penaltyCost = 10;
		double[] meanDemand = { 20, 40, 60, 40 };
		double holdingCost = 1;
		int maxOrderQuantity = 150;
		boolean isForDrawGy = true;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new NormalDist(meanDemand[i], 0.25*meanDemand[i])) // can be changed to other distributions
				.toArray(Distribution[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		
//		int T = 7;
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
		System.out.println("final optimal expected value is: " + recursion.getExpectedValue(initialState));
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");

		/*******************************************************************
		 * Simulating sdp results
		 */
		int sampleNum = 10000;
		Simulation simuation = new Simulation(distributions, sampleNum, recursion);
		simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001; 
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);

		/*******************************************************************
		 * Drawing
		 */
		int minInventorys = 0;
		int maxInventorys = 200; // for drawing pictures
		int xLength = maxInventorys - minInventorys + 1;
		double[][] xQ = new double[xLength][2];
		int index = 0;
		for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
			period = 1;
			xQ[index][0] = initialInventory;
			recursion.getExpectedValue(new State(period, initialInventory));
			xQ[index][1] = recursion.getAction(new State(period, initialInventory));
			index++;
		}
		//Drawing.drawXQ(xQ);

		// since comupteIfAbsent, we need initializing a new class to draw Gy; if not,
		// java would not compute sdp again
		// must redefine stateTransition function and immediate Function;
		StateTransitionFunction<State, Double, Double, State> stateTransition2 = (state, action, randomDemand) -> {
			double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
					: state.getIniInventory() + action - randomDemand;
			nextInventory = nextInventory > maxInventory ? maxInventory : nextInventory;
			nextInventory = nextInventory < minInventory ? minInventory : nextInventory;
			return new State(state.getPeriod() + 1, nextInventory);
		};

		ImmediateValueFunction<State, Double, Double, Double> immediateValue2 = (state, action, randomDemand) -> {
			double fixedCost = 0, variableCost = 0, inventoryLevel = 0, holdingCosts = 0, penaltyCosts = 0;
			if (isForDrawGy == true && state.getPeriod() == 1) {
				fixedCost = 0;
				variableCost = variOrderingCost * state.getIniInventory();
				inventoryLevel = state.getIniInventory() - randomDemand;
			} else {
				fixedCost = action > 0 ? fixedOrderingCost : 0;
				variableCost = variOrderingCost * action;
				inventoryLevel = state.getIniInventory() + action - randomDemand;
			}
			holdingCosts = holdingCost * Math.max(inventoryLevel, 0);
			penaltyCosts = penaltyCost * Math.max(-inventoryLevel, 0);
			double totalCosts = fixedCost + variableCost + holdingCosts + penaltyCosts;
			return totalCosts;
		};

		Recursion recursion2 = new Recursion(OptDirection.MIN, pmf, getFeasibleAction, stateTransition2,
				immediateValue2);
		double[][] yG = new double[xLength][2];
		index = 0;
		for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
			yG[index][0] = initialInventory;
			yG[index][1] = recursion2.getExpectedValue(new State(period, initialInventory));
			index++;
		}
		CheckKConvexity.check(yG, fixedOrderingCost);
		Drawing.drawSimpleG(yG);
		Drawing.drawGAndsS(yG, fixedOrderingCost);
	}
}
