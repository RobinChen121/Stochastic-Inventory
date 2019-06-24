package cash.strongconstraint;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import cash.strongconstraint.FindsCS.FindCCrieria;
import milp.MipCashConstraint;
import sdp.cash.CashRecursion;
import sdp.cash.CashRecursion.OptDirection;
import sdp.cash.CashSimulation;
import sdp.inventory.CheckKConvexity;
import sdp.inventory.GetPmf;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.cash.CashState;
import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date 2018, March 3th, 6:31:10 pm
 * @Description stochastic lot sizing problem with strong cash balance
 *              constraint, provide a (s, C, S) policy
 *
 * a numerical case:
 * double[] meanDemand = {41.8, 6.6, 2, 21.8};
 * double iniCash = 15;
 * double iniInventory = 0;
 * double fixOrderCost = 10;
 * double variCost = 1;
 * double price = 8;
 * double salvageValue = 0.5;
 *
 * there are states: [3, 14, 16, 6], [3, 14, 23, 0]
 */

public class CashConstraint {
	
	// d=[8, 10, 10], iniCash=20, K=10; price=5, v=1; h = 1
	public static void main(String[] args) {
		double[] meanDemand = {15, 15, 15, 15, 15, 15, 15, 15};
		double iniInventory = 0;
		int iniCash = 50; // integer or not does not affect too much
		double fixOrderCost = 15;
		double variCost = 1;
		double holdingCost = 1.5;
		double price = 4;
		double salvageValue = 0;
		FindCCrieria criteria = FindCCrieria.XRELATE; // criteria for finding C			
		double minCashRequired = 0; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;
		
		double discountFactor = 1;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				//.mapToObj(i -> new NormalDist(meanDemand[i], 0.25 * meanDemand[i])) // can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i]))
				.toArray(Distribution[]::new);

//		double[] values1 = {6, 7};
//		double[] probs1 = {0.95, 0.05};
//		double[] values2 = {6, 7};
//		double[] probs2 = {0.95, 0.05};
//		DiscreteDistribution[] distributions = new DiscreteDistribution[T];
//		for (int i = 0; i < T; i++) {
//			if (i % 2 == 0) 
//				distributions[i] = new DiscreteDistribution(values1, probs1, values1.length);
//			else
//				distributions[i] = new DiscreteDistribution(values2, probs2, values2.length);			
//		}	
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
			

		// feasible actions
		Function<CashState, double[]> getFeasibleAction = s -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (s.getIniCash() - minCashRequired - fixOrderCost) / variCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue - fixedCost - variableCost - holdCosts;
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
			nextCash = Math.round(nextCash * 1) / 1; 
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
		double finalValue = iniCash + recursion.getExpectedValue(initialState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		/*******************************************************************
		 * Simulating sdp results
		 */
		int sampleNum = 10000;
		
		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor, 
				fixOrderCost, price, variCost, holdingCost, salvageValue);
		double simFinalValue = simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001; 
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);
		
		/*******************************************************************
		 * Find (s, C1, C2 S) by SDP and simulate
		 */
		System.out.println("");
		double[][] optTable = recursion.getOptTable();
		FindsCS findsCS = new FindsCS(iniCash, distributions, fixOrderCost, price, variCost, holdingCost, salvageValue);
		double[][] optsC12S = findsCS.getsC12S(optTable, minCashRequired, criteria);
		Map<State, Double> cacheC1Values = new TreeMap<>();
		Map<State, Double> cacheC2Values = new TreeMap<>();
		cacheC1Values = findsCS.cacheC1Values;
		cacheC2Values = findsCS.cacheC2Values;
		double simsCSFinalValue = simuation.simulatesCS(initialState, optsC12S, cacheC1Values, cacheC2Values,
				minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
		double gap1 = (finalValue -simsCSFinalValue)/finalValue;
		double gap2 = (simFinalValue -simsCSFinalValue)/simFinalValue;
		System.out.printf("Optimality gap for (s, C1, C2 S) is: %.2f%% or %.2f%%\n", gap1 * 100, gap2 * 100);

		/*******************************************************************
		 * Find (s, C1, S) by SDP and simulate
		 */
		System.out.println("");
		System.out.println("************************************************");
		double[][] optsCS = findsCS.getsCS(optTable, minCashRequired, criteria);
		cacheC1Values = findsCS.cacheC1Values;
		simsCSFinalValue = simuation.simulatesCS(initialState, optsCS, cacheC1Values,
				minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
		double gap21 = (finalValue -simsCSFinalValue)/finalValue;
		double gap22 = (simFinalValue -simsCSFinalValue)/simFinalValue;
		System.out.printf("Optimality gap for (s, C1, S) is: %.2f%% or %.2f%%\n", gap21 * 100, gap22 * 100);
		
		/*******************************************************************
		 * Find (s, meanC, S) by SDP and simulate
		 */
		System.out.println("");
		System.out.println("************************************************");
		optsCS = findsCS.getsCS(optTable, minCashRequired, FindCCrieria.AVG);
		simsCSFinalValue = simuation.simulatesMeanCS(initialState, optsCS, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
		double gap31 = (finalValue -simsCSFinalValue)/finalValue;
		double gap32 = (simFinalValue -simsCSFinalValue)/simFinalValue;
		System.out.printf("Optimality gap for (s, meanC, S) is: %.2f%% or %.2f%%\n", gap31 * 100, gap32 * 100);

//		/*******************************************************************
//		 * Check (s, C1, C2, S) policy, 
//		 * sometimes not always hold, because in certain period 
//		 * for some state C is 12, and 13 in other state, 
//		 * we use heuristic step by choosing maximum one
//		 */		
// 		findsCS.checksC12S(optsC12S, optTable, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
// 		System.out.printf(
//				"\n*******************************************************************\n");
 		
 		/*******************************************************************
		 * Find (s, C, S) by MIP and simulate
		 */
		System.out.println("************************************************");
 		MipCashConstraint mipHeuristic = new MipCashConstraint(iniInventory, iniCash, fixOrderCost, variCost, holdingCost, price, salvageValue, distributions);
 		double[][] sCS = mipHeuristic.findsCS(); 
 		cacheC1Values = mipHeuristic.cacheC1Values;
 		double simsCSMIPValue = simuation.simulatesCS(initialState, sCS, cacheC1Values, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
		gap1 = (finalValue - simsCSMIPValue)/finalValue;
		gap2 = (simFinalValue - simsCSMIPValue)/simFinalValue;	
		System.out.printf("Optimality gap is: %.2f%% or %.2f%%\n", gap1 * 100, gap2 * 100);
		
 		
 		/*******************************************************************
		 * Check K-convexity
		 */	
// 		int minInventorys = 0;
//		int maxInventorys = 100; 
//		int xLength = maxInventorys - minInventorys + 1;
// 		double[][] yG = new double[xLength][2];
//		int index = 0;
//		for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
//			yG[index][0] = initialInventory;
//			yG[index][1] = -recursion.getExpectedValue(new CashState(period, initialInventory, iniCash));
//			index++;
//		}
//		int capacity = (int) Math.max(0, (iniCash - minCashRequired - fixOrderCost) / variCost);
// 		CheckKConvexity.checkCK(yG, fixOrderCost, capacity);
	}
	
	
}
