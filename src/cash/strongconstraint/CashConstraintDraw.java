package cash.strongconstraint;

import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.cash.CashRecursion;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.Drawing;
import sdp.inventory.GetPmf;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;


/**
 * @author chen zhen
 * @version 创建时间：2018年4月15日 上午10:32:46
 * @Description: drawing pictures, see K convexity
 */
public class CashConstraintDraw {

	public static void main(String[] args) {
		double[] meanDemand = { 8, 10, 10 };

		double iniCash = 20;
		double fixOrderCost = 10;
		double variCost = 1;
		double holdingCost = 1;
		double price = 5;
		double minCashRequired = 0; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 2000; // maximum ordering quantity when having enough cash
		boolean isForDrawGy = true;

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other distributions
				.toArray(PoissonDist[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		// feasible actions
		Function<CashState, double[]> getFeasibleAction = s -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (s.getIniCash() - -minCashRequired - fixOrderCost) / variCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = 0;
			double fixedCost = 0;
			double variableCost = 0;
			double inventoryLevel = 0;
			if (isForDrawGy == true && state.getPeriod() == 1) {
				revenue = price * Math.min(state.getIniInventory(), randomDemand);
				fixedCost = 0;
				variableCost = variCost * state.getIniInventory();
				inventoryLevel = state.getIniInventory() - randomDemand;
			} else {
				revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
				fixedCost = action > 0 ? fixOrderCost : 0;
				variableCost = variCost * action;
				inventoryLevel = state.getIniInventory() + action - randomDemand;
			}
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue - fixedCost - variableCost - holdCosts;
			return cashIncrement;
		};

		// state transition function
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition = (state, action,
				randomDemand) -> {
			double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
					: state.getIniInventory() + action - randomDemand;
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double holdCosts = holdingCost * Math.max(nextInventory, 0);
			double nextCash = state.getIniCash() + revenue - fixedCost - variableCost - holdCosts;
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			// cash is integer or not
			// nextCash = Math.round(nextCash * 100) / 100.00;
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};

		/*******************************************************************
		 * Solve
		 */
		CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition,
				immediateValue);
		int period = 1;
		double iniInventory = 0;
		CashState initialState = new CashState(period, iniInventory, iniCash);
		long currTime = System.currentTimeMillis();
		recursion.setTreeMapCacheAction();
		double finalCash = iniCash + recursion.getExpectedValue(initialState);
		System.out.println("final optimal cash is: " + finalCash);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");

		/*******************************************************************
		 * Simulating sdp results
		 */
		int sampleNum = 10000;
		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion);
		simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001;
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);

		/*******************************************************************
		 * Drawing x Q
		 */
		int minInventorys = 0;
		int maxInventorys = 50; // for drawing pictures
		int xLength = maxInventorys - minInventorys + 1;
		double[][] xQ = new double[xLength][2];
		int index = 0;
		for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
			period = 1;
			xQ[index][0] = initialInventory;
			recursion.getExpectedValue(new CashState(period, initialInventory, iniCash));
			xQ[index][1] = recursion.getAction(new CashState(period, initialInventory, iniCash));
			index++;
		}
		Drawing drawing = new Drawing();
		drawing.drawXQ(xQ);

		/*******************************************************************
		 * Drawing y G since comupteIfAbsent, we need initializing a new class to draw
		 * Gy; if not, java would not compute sdp again, we must redefine
		 * stateTransition function and immediate Function;
		 */
		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue2 = (state, action, randomDemand) -> {
			double revenue = 0;
			double fixedCost = 0;
			double variableCost = 0;
			double inventoryLevel = 0;
			if (isForDrawGy == true && state.getPeriod() == 1) {
				revenue = price * Math.min(state.getIniInventory(), randomDemand);
				fixedCost = 0;
				variableCost = variCost * state.getIniInventory();
				inventoryLevel = state.getIniInventory() - randomDemand;
			} else {
				revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
				fixedCost = action > 0 ? fixOrderCost : 0;
				variableCost = variCost * action;
				inventoryLevel = state.getIniInventory() + action - randomDemand;
			}
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue - fixedCost - variableCost - holdCosts;
			return cashIncrement;
		};

		// state transition function
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition2 = (state, action,
				randomDemand) -> {
			double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
					: state.getIniInventory() + action - randomDemand;
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double holdCosts = holdingCost * Math.max(nextInventory, 0);
			double nextCash = state.getIniCash() + revenue - fixedCost - variableCost - holdCosts;
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};

		CashRecursion recursion2 = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition2,
				immediateValue2);
		double[][] yG = new double[xLength][2];
		index = 0;
		for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
			yG[index][0] = initialInventory;
			yG[index][1] = -recursion2.getExpectedValue(new CashState(period, initialInventory, iniCash));
			index++;
		}
		drawing.drawSimpleG(yG);
		drawing.drawGAndsS(yG, fixOrderCost);
	}

}
