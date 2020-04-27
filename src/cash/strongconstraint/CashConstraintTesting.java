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
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.CheckKConvexity;
import sdp.inventory.Drawing;
import sdp.inventory.GetPmf;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToCsv;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author chen zhen
 * @version 2018, March 3rd, 6:31:10 pm
 * @Description: strong cash constraints testing
 * 
 * there are some states, same initial inventory, different initial cash,
 * low initial cash order while high initial cash not ordering.
 * 
 * Check CK-convexity, GB not CK-convex, but GA seems to be CK-convex
 * 
 * MIP model performs better when planning horizon length is longer
 */

public class CashConstraintTesting {

	// average computation time for 10 periods is 500s, 9 periods is 150s, or 305s,
	// or 400s
	public static void main(String[] args) {
		String headString = "K, v, h, I0, price, salvageValue, B0, DemandPatt, OptValue, Time(sec), simValue, "
				+ "totalStates, "
				+ "simsC1SValue, gap1, gap2, "
				+ "simsMeanCSValue, gap1, gap2, "
				+ "firstQ, capacity, mipValue, gap11, gap22, time2";
		WriteToCsv.writeToFile("./" + "test_results.csv", headString);

		double[][] iniMeanDemands = { { 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 },
				{ 21.15, 18.9, 17.7, 16.5, 15.15, 13.95, 12.75, 11.55, 10.35, 9.15}, { 6.6, 9.3, 11.1, 12.9, 16.8, 21.6, 24, 26.4, 31.5, 33.9 },
				{ 12.1,10,7.9,7,7.9,10,12.1,13,12.1,10 }, {15.7,10,4.3,2,4.3,10,15.7,18,15.7,10 },
				{ 41.8, 6.6, 2, 21.8, 44.8, 9.6, 2.6, 17, 30, 35.4 }, { 4.08, 12.16, 37.36, 21.44, 39.12, 35.68, 19.84, 22.48, 29.04, 12.4 },
				{ 4.7, 8.1, 23.6, 39.4, 16.4, 28.7, 50.8, 39.1, 75.4, 69.4 }, { 4.4, 11.6, 26.4, 14.4, 14.6, 19.8, 7.4, 18.3, 20.4, 11.4 },
				{ 4.9, 18.8, 6.4, 27.9, 45.3, 22.4, 22.3, 51.7, 29.1, 54.7 } };

		double[] K = {10, 15, 20};
		double[] v = {1};
		double[] B0 = { 3, 5, 7}; // ini cash can order 4 or 6 items
		double[] p = { 5, 6, 7};  // margin is 4, 5, 6
		double[] h = {0};
		double salvageValue = 0;	
		
		FindCCrieria criteria = FindCCrieria.XRELATE;
		double truncationQuantile = 0.999;
		int stepSize = 1;
		double overheadCost = 0; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 150; // maximum ordering quantity when having enough cash
		double minInventoryState = 0;
		double maxInventoryState = 200;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 1500;
		double discountFactor = 1; // generally 1 for the cash constrained problem

		/*******************************************************************
		 * set demands length, for testing 
		 */
		int newLength = iniMeanDemands[0].length;
		double[][] meanDemands = new double[iniMeanDemands.length][newLength];
		for (int i = 0; i < iniMeanDemands.length; i++)
			for (int j = 0; j < newLength; j++) {
				meanDemands[i][j] = iniMeanDemands[i][j];
			}
		
		 
		for (int idemand = 0; idemand < meanDemands.length; idemand++)
			for (int iK = 0; iK < K.length; iK++)
				for (int iv = 0; iv < v.length; iv++)
					for (int ip = 0; ip < p.length; ip++)
						for (int ih = 0; ih < h.length; ih++)
							for (int iB = 0; iB < B0.length; iB++) {
								double[] meanDemand = meanDemands[idemand];
								double fixOrderCost = K[iK];
								double variCost = v[iv];
								double price = p[ip];
								double iniCash = fixOrderCost + overheadCost + variCost * B0[iB];
								double holdingCost = h[ih];
								double iniInventory = 0;
								
								
								// get demand possibilities for each period
								int T = meanDemand.length;
								Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
										//.mapToObj(i -> new NormalDist(meanDemand[i], 0.25 * meanDemand[i])) // can be changed to other distributions
										.mapToObj(i -> new PoissonDist(meanDemand[i]))
										.toArray(Distribution[]::new);
								double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

								// feasible actions
								Function<CashState, double[]> getFeasibleAction = s -> {
									double maxQ = (int) Math.min(maxOrderQuantity,
											Math.max(0, (s.getIniCash() - overheadCost - fixOrderCost) / variCost));
									return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
								};

								// immediate value
								ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state,
										action, randomDemand) -> {
									double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
									double fixedCost = action > 0 ? fixOrderCost : 0;
									double variableCost = variCost * action;
									double inventoryLevel = state.getIniInventory() + action - randomDemand;
									double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
									double cashIncrement = revenue - fixedCost - variableCost - holdCosts - overheadCost;
									double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
									cashIncrement += salValue;
									return cashIncrement;
								};

								// state transition function
								StateTransitionFunction<CashState, Double, Double, CashState> stateTransition = (state,
										action, randomDemand) -> {
									double nextInventory = Math.max(0, state.getIniInventory() + action - randomDemand);
									double nextCash = state.getIniCash()  + immediateValue.apply(state, action, randomDemand);;
									nextCash = nextCash > maxCashState ? maxCashState : nextCash;
									nextCash = nextCash < minCashState ? minCashState : nextCash;
									nextInventory = nextInventory > maxInventoryState ? maxInventoryState
											: nextInventory;
									nextInventory = nextInventory < minInventoryState ? minInventoryState
											: nextInventory;
									// cash is integer or not
									nextCash = Math.round(nextCash * 1) / 1;
									return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
								};

								/*******************************************************************
								 * Solve
								 */
								CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction,
										stateTransition, immediateValue, discountFactor);
								int period = 1;

								CashState initialState = new CashState(period, iniInventory, iniCash);
								long currTime = System.currentTimeMillis();
								recursion.setTreeMapCacheAction();
								double finalValue = recursion.getExpectedValue(initialState) + iniCash;
								System.out.println("final optimal expected cash is: " + finalValue);
								double firstQ = recursion.getAction(initialState);
								System.out.println("optimal order quantity in the first priod is : " + firstQ);
								double time = (System.currentTimeMillis() - currTime) / 1000.0;
								System.out.println("running time is " + time + "s");

								/*******************************************************************
								 * Simulating sdp results
								 */
								int sampleNum = 100000;
								CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion,
										discountFactor, fixOrderCost, price, variCost, holdingCost, salvageValue);
								double simFinalValue = simuation.simulateSDPGivenSamplNum(initialState);
								
								/*******************************************************************
								 * Find (s, C1, C2, S) and simulate
								 */
								System.out.println("");
								System.out.println("************************************************");
								double[][] optTable = recursion.getOptTable();
								FindsCS findsCS = new FindsCS(iniCash, distributions, fixOrderCost, price, variCost, holdingCost, salvageValue);
//								double[][] optsCS = findsCS.getsC12S(optTable, overheadCost, criteria);
								Map<State, Double> cacheC1Values = new TreeMap<>();
//								Map<State, Double> cacheC2Values = new TreeMap<>();
//								cacheC1Values = findsCS.cacheC1Values;
//								cacheC2Values = findsCS.cacheC2Values;
//								double simsC1C2SFinalValue = simuation.simulatesCS(initialState, optsCS, cacheC1Values, cacheC2Values,
//										overheadCost, maxOrderQuantity, fixOrderCost, variCost);
//								double gapsC1C2S1 = (finalValue - simsC1C2SFinalValue) / finalValue;
//								double gapsC1C2S2 = (simFinalValue - simsC1C2SFinalValue) / simFinalValue;
//								System.out.printf("Optimality gap for (s, C1, C2, S) is: %.2f%% or %.2f%%\n", gapsC1C2S1 * 100, gapsC1C2S2 * 100);
																
								
								/*******************************************************************
								 * Find (s, C1, S) and simulate
								 */
								System.out.println("");
								System.out.println("************************************************");
								double[][] optsCS = findsCS.getsCS(optTable, overheadCost, criteria);
								cacheC1Values = findsCS.cacheC1Values;
								double simsC1SFinalValue = simuation.simulatesCS(initialState, optsCS, cacheC1Values,
										overheadCost, maxOrderQuantity, fixOrderCost, variCost);
								double gapsC1S1 = (finalValue - simsC1SFinalValue) / finalValue;
								double gapsC1S2 = (simFinalValue - simsC1SFinalValue) / simFinalValue;
								System.out.printf("Optimality gap for (s, C1, S) is: %.2f%% or %.2f%%\n", gapsC1S1 * 100, gapsC1S2 * 100);
								
								/*******************************************************************
								 * Find (s, meanC, S) by SDP and simulate
								 */
								System.out.println("");
								System.out.println("************************************************");							
								optsCS = findsCS.getsCS(optTable, overheadCost, FindCCrieria.AVG);
								double simsMeanCSFinalValue = simuation.simulatesMeanCS(initialState, optsCS, overheadCost, maxOrderQuantity, fixOrderCost, variCost);
								double gapsMeanCS1 = (finalValue - simsMeanCSFinalValue)/finalValue;
								double gapsMeanCS2 = (simFinalValue - simsMeanCSFinalValue)/simFinalValue;
								System.out.printf("Optimality gap for (s, meanC, S) is: %.2f%% or %.2f%%\n", gapsMeanCS1 * 100, gapsMeanCS2 * 100);
								

//								/*******************************************************************
//								 * Check (s, C, S) policy, sometimes not always hold
//								 */
//								int nonOptCount = findsCS.checksCS(optsCS, optTable, overheadCost, maxOrderQuantity,
//										fixOrderCost, variCost);
								
						 		/*******************************************************************
								 * Find (s, C, S) by MIP and simulate
								 */
								double simsCSMIPValue = 0;
								double gap11 = 0;
								double gap22 = 0;
								double time2 = 0;
								currTime = System.currentTimeMillis();
						 		MipCashConstraint mipHeuristic = new MipCashConstraint(iniInventory, iniCash, fixOrderCost, variCost, holdingCost, price, salvageValue, distributions, overheadCost);
						 		double[][] sCS = mipHeuristic.findsCS(); 					 		
						 		time2 = (System.currentTimeMillis() - currTime) / 1000.0;
								System.out.println("running time is " + time2 + "s");
						 		cacheC1Values = mipHeuristic.cacheC1Values;
						 		simsCSMIPValue = simuation.simulatesCS(initialState, sCS, cacheC1Values, overheadCost, maxOrderQuantity, fixOrderCost, variCost);
								gap11 = (finalValue - simsCSMIPValue)/finalValue;
								gap22 = (simFinalValue - simsCSMIPValue)/simFinalValue;	
								System.out.printf("Optimality gap by Mip is: %.2f%% or %.2f%%\n", gap11 * 100, gap22 * 100);
								
								System.out.printf(
										"\n*******************************************************************\n");
								
								/*******************************************************************
								 * Check K-convexity of GA or GB
								 */
								// immediate value
//								boolean isForDrawGy = true;
//								ImmediateValueFunction<CashState, Double, Double, Double> immediateValue2 = (state, action, randomDemand) -> {
//									double revenue = 0;
//									double fixedCost = 0;
//									double variableCost = 0;
//									double inventoryLevel = 0;
//									if (isForDrawGy == true && state.getPeriod() == 1) {
//										revenue = price * Math.min(state.getIniInventory(), randomDemand);
//										fixedCost = 0;
//										variableCost = variCost * state.getIniInventory();
//										inventoryLevel = state.getIniInventory() - randomDemand;
//									} else {
//										revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
//										fixedCost = action > 0 ? fixOrderCost : 0;
//										variableCost = variCost * action;
//										inventoryLevel = state.getIniInventory() + action - randomDemand;
//									}
//									double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
//									double cashIncrement = revenue - fixedCost - variableCost - holdCosts - overheadCost;
//									double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
//									cashIncrement += salValue;
//									return cashIncrement;
//								};
//								StateTransitionFunction<CashState, Double, Double, CashState> stateTransition2 = (state, action,
//										randomDemand) -> {
//									double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
//												: state.getIniInventory() + action - randomDemand;
//									double nextCash = state.getIniCash() + immediateValue2.apply(state, action, randomDemand);
//									nextCash = nextCash > maxCashState ? maxCashState : nextCash;
//									nextCash = nextCash < minCashState ? minCashState : nextCash;
//									nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
//									nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
//									return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
//								};
//								CashRecursion recursion2 = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition2,
//										immediateValue2, discountFactor);								
//								
//								int minInventorys = 0;
//								int maxInventorys = 50; // for drawing pictures
//								int xLength = maxInventorys - minInventorys + 1;
//								int index = 0;
//								int capacity = (int) Math.max(0, (iniCash - overheadCost - fixOrderCost) / variCost);
//								double[][] yG2 = new double[xLength][2];
//								index = 0;
////								for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
////									yG2[index][0] = initialInventory;
////									yG2[index][1] = -recursion2.getExpectedValue(new CashState(period, initialInventory, iniCash));
////									index++;
////								}		
////						 		
////								Drawing drawing = new Drawing();
////								drawing.drawSimpleG(yG2, iniCash, "K transfered in cash GB"); // GB
////								System.out.println();
////								System.out.println("CK convex of GB:");
////						 		CheckKConvexity.checkCK(yG2, fixOrderCost, capacity); // check the CK convex of GB
////						 		System.out.println();
//						 		
//						 	    // state transition function 3
//						 		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition3 = (state, action,
//						 				randomDemand) -> {
//						 		double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
//						 					: state.getIniInventory() + action - randomDemand;
//						 			double nextCash = state.getIniCash() + immediateValue2.apply(state, action, randomDemand);
//						 			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
//						 			nextCash = nextCash < minCashState ? minCashState : nextCash;
//						 			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
//						 			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
//						 			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
//						 		};
//						 		CashRecursion recursion3 = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition3,
//										immediateValue2, discountFactor);
//								double[][] yG3 = new double[xLength][2];
//								index = 0;
//								for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
//									yG3[index][0] = initialInventory;
//									yG3[index][1] = -recursion3.getExpectedValue(new CashState(period, initialInventory, iniCash));
//									index++;
//								}
//								
								String convexity = "";
//								double[][] yG = new double[xLength][2];
//								index = 0;
//								for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
//									yG[index][0] = initialInventory;
//									yG[index][1] = -recursion2
//											.getExpectedValue(new CashState(period, initialInventory, iniCash));
//									index++;
//								}
//								convexity = CheckKConvexity.checkCK(yG3, fixOrderCost, capacity);
//								System.out.printf( 
//										"\n*******************************************************************\n");
								
								long totalStates = optTable.length;
								String out = fixOrderCost + ",\t" + variCost + ",\t" + holdingCost + ",\t"
										+ iniInventory + ",\t" + price + ",\t" + salvageValue + ",\t" + iniCash + ",\t" + (idemand + 1) + ",\t"
										+ finalValue + ",\t" + time + ",\t" + simFinalValue + ",\t"
										+ totalStates + ",\t"+ simsC1SFinalValue +",\t" + gapsC1S1 + ",\t" + gapsC1S2 + ",\t" 
										+ simsMeanCSFinalValue +",\t" + gapsMeanCS1 + ",\t" + gapsMeanCS2 + ",\t" 
										+ firstQ + ",\t" + B0[iB] + ",\t"+ 
										simsCSMIPValue + ",\t" + gap11 + ",\t" + gap22+ ",\t" + time2;

								WriteToCsv.writeToFile("./" + "test_results.csv", out);
							}
	}

}
