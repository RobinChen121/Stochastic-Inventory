package cash.overdraft;

import java.util.Arrays;
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
 *               and near optimal.
 *               otherwise, optimal policy is complex: may have another S level, (s1, S1, C1, C2, s2, S2) policy.
 *               maybe state(x, R) is better to characterize the optimal policy
 *               
 *               
 */

public class CashOverdraft {

	public static void main(String[] args) {
		double[] meanDemand = {10, 20, 10};
		
		double[] overheadCost = {50, 50, 50};
		double fixOrderCost = 0;
		double variCost = 1;
		double holdingCost = 0;
		double price = 5;
		double salvageValue = 0.5;

		double iniCash = 0;
		double r0 = 0.01;
		double r1 = 0;
		double r2 = 0.1;
		double r3 = 2; // penalty interest rate for overdraft exceeding the limit
		double limit = 80; // overdraft limit
		double interestFreeQuantity = 25;
		double maxOrderQuantity = 100; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 100;
		double minCashState = -500;
		double maxCashState = 1000;
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
					Math.max(0, (s.getIniCash() -minCashRequired - fixOrderCost) / variCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		// when to pay interest is important
		// hold cost is zero
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashBalanceBefore = state.getIniCash() - fixedCost - variableCost - holdCosts;// whether plus revenue in this time point
			double interest = interestRate * Math.max(-cashBalanceBefore, 0);
			double cashBalanceAfter = cashBalanceBefore - interest + revenue;
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
		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor);
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
//		System.out.println("(s, S) are: ");
//		System.out.println(Arrays.deepToString(optsS));
		double simsSFinalValue = simuation.simulatesSOD(initialState, optsS, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
		System.out.printf("Optimality gap for policy (s, S) is: %.2f%%\n", (finalCash -simsSFinalValue)/finalCash*100);
		
		
		/*******************************************************************
		 * Find (s, C, S1, S2) and simulate
		 * 
		 */
		double[][] optsCS1S2 = findsS.getsCS1S2(optTable);
		double simsSFinalValue2 = simuation.simulatesSOD2(initialState, optsCS1S2);
		System.out.printf("Optimality gap for policy (s, C, S1, S2) is: %.2f%%\n", (finalCash -simsSFinalValue2)/finalCash*100);
		
		
		
		
	}

}
