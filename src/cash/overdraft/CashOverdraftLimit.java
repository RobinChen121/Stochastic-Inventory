package cash.overdraft;

import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.cash.CashRecursion;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author chen zhen
 * @version 2018, April 10th, 9:38:22 pm
 * @Description: consider business overdraft with limit in a cash stochastic lot
 *               sizing problem
 */
public class CashOverdraftLimit {

	public static void main(String[] args) {
		double[] meanDemand = {10, 15, 10};
		
		double[] overheadCost = {50, 50, 50};

		double fixOrderCost = 0;
		double variCost = 1;
		double holdingCost = 0;
		double price = 10;
		double salvageValue = 0.5;

		double iniCash = 0;
		double interestRate = 0.1; // about 10 times higher than the deposite rate
		double depositeRate = 0.01;
		double minCashRequired = -100; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 100; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 100;
		double minCashState = -200; // 
		double maxCashState = 800;
		double discountFactor = 1;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other distributions
				.toArray(PoissonDist[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		// feasible actions
		Function<CashState, double[]> getFeasibleAction = s -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (s.getIniCash() -overheadCost[s.getPeriod()-1] - minCashRequired - fixOrderCost) / variCost));
			maxQ = maxOrderQuantity;
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashBalanceBeforeRevenue = state.getIniCash() - fixedCost - variableCost - holdCosts - overheadCost[state.getPeriod()-1];
			// related with when to pay interest
			double interest = interestRate * Math.max(-cashBalanceBeforeRevenue, 0); 
			double deposite = depositeRate * Math.max(cashBalanceBeforeRevenue, 0);
			double cashBalanceAfter = cashBalanceBeforeRevenue - interest + deposite + revenue;
			double cashIncrement = cashBalanceAfter - state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
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
			//cash is integer or not
			nextCash = Math.round(nextCash * 10) / 10;
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};

		/*******************************************************************
		 * Solve
		 */
		CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition,
				immediateValue, discountFactor);
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
		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor);
		simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001; 
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);

		/*******************************************************************
		 * Find (s, C, S) and simulate
		 */
		System.out.println("");
		double[][] optTable = recursion.getOptTable();
		FindsSOverDraft findsCS = new FindsSOverDraft(T, iniCash);
//		double[][] optsCS = findsCS.getsCS(optTable);		
//		double simsCSFinalValue = simuation.simulatesCSDraft(initialState, optsCS, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
//		System.out.printf("Optimality gap is: %.2f%%\n", (finalCash -simsCSFinalValue)/finalCash*100);
	}
}
