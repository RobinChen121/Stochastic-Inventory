package cash.strongconstraint;

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
import sdp.write.WriteToCsv;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/** 
 * @author chen zhen 
 * @version 2018, March 3rd, 6:31:10 pm 
 * @Description: strong cash constraints testing
 */

public class CashConstraintTesting {

	// average computation time for 10 periods is 500s, 9 periods is 150s, or 305s, or 400s
	public static void main(String[] args) {
		String headString = "K, v, h, I0, price, B0, DemandPatt, OptValue, Time(sec), simValue, simsCSValue, nonOptStatesCount, gap1, gap2";
		WriteToCsv.writeToFile("./" + "test_results.csv", headString);

		double[][] meanDemands = {{15, 15, 15, 15, 15, 15, 15, 15},
				{21.15, 18.9, 17.7, 16.5, 15.15, 13.95, 12.75, 11.55},
				{6.6, 9.3, 11.1, 12.9, 16.8, 21.6, 24, 26.4},
				{9, 13, 20, 16, 10, 16, 22, 15},
				{22, 18, 11, 16, 22, 12, 14, 21},
				{41.8, 6.6, 0.4, 21.8, 44.8, 9.6, 2.6, 17},							
				{4.08, 12.16, 37.36, 21.44, 39.12, 35.68, 19.84, 22.48},
				{4.7, 8.1, 23.6, 39.4, 16.4, 28.7, 50.8, 39.1},
				{4.4, 11.6, 26.4, 14.4, 14.6, 19.8, 7.4, 18.3},
				{4.9, 18.8, 6.4, 27.9, 45.3, 22.4, 22.3, 51.7}};

		double[] K = {10, 20}; 
		double[] v = {1, 2};
		double[] B0 = {1, 2}; // change to K + v * 5
		double[] p = {4, 8};

		double truncationQuantile = 0.9999;  
		int stepSize = 1;
		double holdingCost = 1;
		double minCashRequired = 0; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash		
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -100;  // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;			

		for (int idemand = 0; idemand < meanDemands.length; idemand++) 
			for (int iK = 0; iK < K.length; iK++) 
				for (int iv = 0; iv < v.length; iv++) 
					for (int ip = 0; ip < p.length; ip++)
						for (int iB = 0; iB < B0.length; iB++) { 
							double[] meanDemand = meanDemands[idemand];						
							double fixOrderCost = K[iK];		
							double variCost = v[iv];		
							double price = p[ip];	
							double iniCash = fixOrderCost + variCost * B0[iB] * 5; 
							
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
								double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
								double fixedCost = action > 0 ? fixOrderCost : 0;
								double variableCost = variCost * action;
								double inventoryLevel = state.getIniInventory() + action - randomDemand;
								double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
								double cashIncrement = revenue - fixedCost - variableCost - holdCosts;
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
										nextCash = nextCash > maxCashState ? maxCashState : nextCash;
										nextCash = nextCash < minCashState ? minCashState : nextCash;
										nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
										nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
										// cash is integer or not
										//nextCash = Math.round(nextCash * 100) / 100.00; 
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
									double finalValue = recursion.getExpectedValue(initialState) + iniCash;
									System.out.println("final optimal expected cash is: " + finalValue);
									System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
									double time = (System.currentTimeMillis() - currTime) / 1000;
									System.out.println("running time is " + time + "s");

									/*******************************************************************
									 * Simulating sdp results
									 */
									int sampleNum = 10000;
									CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion);
									double simFinalValue = simuation.simulateSDPGivenSamplNum(initialState);
									
									/*******************************************************************
									 * Find (s, C, S) and simulate
									 */
									System.out.println("");
									double[][] optTable = recursion.getOptTable();
									FindsCS findsCS = new FindsCS(T, iniCash);
									double[][] optsCS = findsCS.getsCS(optTable);
									double simsCSFinalValue = simuation.simulatesCS(initialState, optsCS, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
									double gap1 = (finalValue -simsCSFinalValue)/finalValue;
									double gap2 = (simFinalValue -simsCSFinalValue)/simFinalValue;	
									System.out.printf("Optimality gap is: %.2f%% or %.2f%%\n", gap1 * 100, gap2 * 100);
									
									/*******************************************************************
									 * Check (s, C, S) policy, 
									 * sometimes not always hold, because in certain period 
									 * for some state C is 12, and 13 in other state, 
									 * we use heuristic step by choosing maximum one
									 */		
									int nonOptCount = findsCS.checksBS(optsCS, optTable, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
									System.out.printf("\n*******************************************************************\n");
														
									String out = fixOrderCost + ",\t"+
											variCost + ",\t"+
											holdingCost + ",\t"+
											iniInventory + ",\t"+
											price + ",\t" +
											iniCash + ",\t" +
											(idemand + 1) + ",\t" +
											finalValue + ",\t" +
											time + ",\t" + 
											simFinalValue + ",\t" +
											simsCSFinalValue + ",\t" +
											nonOptCount + ",\t" +
											gap1 + ",\t" +
											gap2;

									WriteToCsv.writeToFile("./" + "test_results.csv", out);
						}
	}

}
