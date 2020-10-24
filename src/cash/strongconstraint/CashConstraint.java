package cash.strongconstraint;

import java.util.ArrayList;
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
import sdp.inventory.Drawing;
import sdp.inventory.GetPmf;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToCsv;
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
 *              constraint, provide a (s, C, S) policy,
 *              when there is no fixed ordering cost, it is a similar base-stock
 *              policy. 
 * when no inventory holding cost, (s, C, S) may not be optimal, because sometimes may be exist C2,
 * or S is not always follows, S may be x during the states
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
		double[] meanDemand = {15, 15, 15, 15};
		//double[] meanDemand = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
		double iniInventory = 0;
		double iniCash = 13;
		double fixOrderCost = 10;
		double variCost = 1;
		double price = 5;
		double depositeRate = 0;
		double salvageValue = 0.5;
		double holdingCost = 0.5;	
		FindCCrieria criteria = FindCCrieria.XRELATE;		
		double overheadCost = 0; // costs like wages or rents which is required to pay in each period
		double overheadRate = 0; // rate from revenue to pay overhead wages
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;
		double penaltyCost = 10000;
		
		double discountFactor = 1;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				//.mapToObj(i -> new NormalDist(meanDemand[i], Math.sqrt(meanDemand[i]))) // can be changed to other distributions
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
					Math.max(0, (s.getIniCash() - overheadCost - fixOrderCost) / variCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double deposite = (state.getIniCash() - fixedCost - variableCost) * (1 + depositeRate);
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = (1 - overheadRate)*revenue + deposite - holdCosts - overheadCost - state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			double endCash = state.getIniCash() + cashIncrement;
			if (endCash < 0) {
				cashIncrement += penaltyCost * endCash;
			}	
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
		 * parameter vales like price, variCost, holdingCost etc.
		 * are only for compute L(y), not very necessary
		 */
		int sampleNum = 10000;
		
		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor, 
				fixOrderCost, price, variCost, holdingCost, salvageValue); // no need to add overheadCost in this class
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
		double[][] optsC12S = findsCS.getsC12S(optTable, overheadCost, criteria);
		Map<State, Double> cacheC1Values = new TreeMap<>();
		Map<State, Double> cacheC2Values = new TreeMap<>();
 		
		
		double simsCSFinalValue = simuation.simulatesCS(initialState, optsC12S, cacheC1Values, cacheC2Values,
				overheadCost, maxOrderQuantity, fixOrderCost, variCost);
		double gap1 = (finalValue -simsCSFinalValue)/finalValue;
		double gap2 = (simFinalValue -simsCSFinalValue)/simFinalValue;
		System.out.printf("Optimality gap for (s, C1, C2 S) is: %.2f%% or %.2f%%\n", gap1 * 100, gap2 * 100);

		/*******************************************************************
		 * Find (s, C1, S) by SDP and simulate
		 */
		System.out.println("");
		System.out.println("************************************************");
		double[][] optsCS = findsCS.getsCS(optTable, overheadCost, criteria);
		cacheC1Values = findsCS.cacheC1Values;
		simsCSFinalValue = simuation.simulatesCS(initialState, optsCS, cacheC1Values,
				overheadCost, maxOrderQuantity, fixOrderCost, variCost);
		double gap21 = (finalValue -simsCSFinalValue)/finalValue;
		double gap22 = (simFinalValue -simsCSFinalValue)/simFinalValue;
		System.out.printf("Optimality gap for (s, C1, S) is: %.2f%% or %.2f%%\n", gap21 * 100, gap22 * 100);
		//double[][] numFrequency = findsCS.getMaxSFrequency(optTable, overheadCost, criteria);
		//System.out.println("most frequent S in each period");
		//System.out.println(Arrays.deepToString(numFrequency));
		
		
		/*******************************************************************
		 * Find (s, meanC, S) by SDP and simulate
		 */
		System.out.println("");
		System.out.println("************************************************");
		optsCS = findsCS.getsCS(optTable, overheadCost, FindCCrieria.AVG);
		simsCSFinalValue = simuation.simulatesMeanCS(initialState, optsCS, overheadCost, maxOrderQuantity, fixOrderCost, variCost);
		double gap31 = (finalValue -simsCSFinalValue)/finalValue;
		double gap32 = (simFinalValue -simsCSFinalValue)/simFinalValue;
		System.out.printf("Optimality gap for (s, meanC, S) is: %.2f%% or %.2f%%\n", gap31 * 100, gap32 * 100);

		/*******************************************************************
		 * Check (s, C1, C2, S) policy, 
		 * sometimes not always hold, because in certain period 
		 * for some state C is 12, and 13 in other state, 
		 * we use heuristic step by choosing maximum one
		 */		
// 		findsCS.checksCS(optsCS, optTable, overheadCost, maxOrderQuantity, fixOrderCost, variCost);
// 		System.out.printf(
//				"\n*******************************************************************\n");
 		
 		/*******************************************************************
		 * Find (s, C, S) by MIP and simulate
		 */
		System.out.println("************************************************");
 		MipCashConstraint mipHeuristic = new MipCashConstraint(iniInventory, iniCash, fixOrderCost, variCost, holdingCost, price, salvageValue, distributions, overheadCost);
 		double[][] sCS = mipHeuristic.findsCSNew(); 
 		Map<State, Double> cacheCValues = new TreeMap<>();
 		cacheCValues = mipHeuristic.cacheC1Values;
 		double simsCSMIPValue = simuation.simulatesCS(initialState, sCS, cacheCValues, overheadCost, maxOrderQuantity, fixOrderCost, variCost);
		gap1 = (finalValue - simsCSMIPValue)/finalValue;
		gap2 = (simFinalValue - simsCSMIPValue)/simFinalValue;	
		System.out.printf("Optimality gap is: %.2f%% or %.2f%%\n", gap1 * 100, gap2 * 100);
		
 		
// 		/*******************************************************************
//		 * Check CK-convexity
//		 */	
//		
//		// immediate value
//		boolean isForDrawGy = true;
//		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue2 = (state, action, randomDemand) -> {
//			double revenue = 0;
//			double fixedCost = 0;
//			double variableCost = 0;
//			double inventoryLevel = 0;
//			if (isForDrawGy == true && state.getPeriod() == 1) {
//				revenue = price * Math.min(state.getIniInventory(), randomDemand);
//				fixedCost = 0;
//				variableCost = variCost * state.getIniInventory();
//				inventoryLevel = state.getIniInventory() - randomDemand;
//			} else {
//				revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
//				fixedCost = action > 0 ? fixOrderCost : 0;
//				variableCost = variCost * action;
//				inventoryLevel = state.getIniInventory() + action - randomDemand;
//			}
//			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
//			double cashIncrement = revenue - fixedCost - variableCost - holdCosts - overheadCost;
//			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
//			cashIncrement += salValue;
//			return cashIncrement;
//		};
//
//		// state transition function 2
//		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition2 = (state, action,
//				randomDemand) -> {
//			double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
//					: state.getIniInventory() + action - randomDemand;
//			double nextCash = state.getIniCash() + immediateValue2.apply(state, action, randomDemand);
//			if (isForDrawGy == true && state.getPeriod() == 1) // this is the only difference with transition function3
//				nextCash -= fixOrderCost;
//			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
//			nextCash = nextCash < minCashState ? minCashState : nextCash;
//			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
//			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
//			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
//		};
//		
//		
//		int minInventorys = 0;
//		int maxInventorys = 50; // for drawing pictures
//		int xLength = maxInventorys - minInventorys + 1;
//		int index = 0;
//		CashRecursion recursion2 = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition2,
//				immediateValue2, discountFactor);
//		double[][] yG2 = new double[xLength][2];
//		index = 0;
//		for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
//			yG2[index][0] = initialInventory;
//			yG2[index][1] = -recursion2.getExpectedValue(new CashState(period, initialInventory, iniCash));
//			index++;
//		}		
// 		
//		int capacity = (int) Math.max(0, (iniCash - overheadCost - fixOrderCost) / variCost);
//		Drawing drawing = new Drawing();
//		drawing.drawSimpleG(yG2, iniCash, "K transfered in cash GB"); // GB
//		System.out.println();
//		System.out.println("CK convex of GB:");
// 		CheckKConvexity.checkCK(yG2, fixOrderCost, capacity); // check the CK convex of GB
// 		System.out.println();
// 		

		
	}			
}
