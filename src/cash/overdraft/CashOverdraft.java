package cash.overdraft;

import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

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
 * @author chen zhen
 * @version 2018, March 15th, 6:31:10 pm
 * @Description: consider overdraft, and lost sales; can't prove K convexity
 *               since interest rate is a disturbing factor; when iniCash is
 *               large enough or small enough, policy (s, S) perform very well
 *               and near optimal otherwise, optimal policy is complex
 */

public class CashOverdraft {

	public static void main(String[] args) {
		double[] meanDemand = { 8, 10, 20};

		double fixOrderCost = 10;
		double variCost = 1;
		double holdingCost = 0;
		double price = 10;
		double salvageValue = 0.5;

		double iniCash = 10;
		double interestRate = 0.8;
		double minCashRequired = -1000; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 100; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 100;
		double minCashState = -500;
		double maxCashState = 1000;
		double discountFactor = 0.95;
		

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other distributions
				.toArray(PoissonDist[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		// feasible actions
		Function<CashState, double[]> getFeasibleAction = s -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (s.getIniCash() -minCashRequired - fixOrderCost) / variCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashBalanceBefore = state.getIniCash() + revenue - fixedCost - variableCost - holdCosts;
			double interest = interestRate * Math.max(-cashBalanceBefore, 0);
			double cashBalanceAfter = cashBalanceBefore - interest;
			double cashIncrement = cashBalanceAfter - state.getIniCash();
			return cashIncrement;
		};

		// state transition function
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition = (state, action,
				randomDemand) -> {
			double nextInventory = Math.max(0, state.getIniInventory() + action - randomDemand);
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double holdCosts = holdingCost * Math.max(nextInventory, 0);
			double nextCash = state.getIniCash() + revenue - fixedCost - variableCost - holdCosts;
			nextCash = nextCash - Math.max(-nextCash, 0) * interestRate;
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
				immediateValue, discountFactor);
		int period = 1; double iniInventory = 0;
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
		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor, 
				fixOrderCost, price, variCost, holdingCost, salvageValue);
		simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001; 
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);

		/*******************************************************************
		 * Find (s, S) and simulate
		 */
		System.out.println("");
		double[][] optTable = recursion.getOptTable();
		FindsSOverDraft findsS = new FindsSOverDraft(T, iniCash);
		double[][] optsS = findsS.getsS(optTable);
		double simsSFinalValue = simuation.simulatesSOD(initialState, optsS, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
		System.out.printf("Optimality gap is: %.2f%%\n", (finalCash -simsSFinalValue)/finalCash*100);
	}

}
