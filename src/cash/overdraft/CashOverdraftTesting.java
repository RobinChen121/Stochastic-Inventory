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

import sdp.write.WriteToCsv;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;


/**
 * @author chen zhen
 * @version 创建时间：2018年3月15日 下午6:31:10
 * @Description: testing overdraft sdp, 
 * if cash is integer, computation time is 35s; if one decimal, time is 312s
 */

public class CashOverdraftTesting {

	public static void main(String[] args) {
		String headString = "K, v, h, I0, pai, B0, rate, DemandPatt, OpValue, Time(sec), simsSValue, gap";
		WriteToCsv.writeToFile("./" + "test_results.csv", headString);

		double[] K = { 10, 15 };
		double[] v = { 1, 2 };
		double[] P = { 5, 10 };
		double[] b = { 0.05, 0.2 };
		int initialCash[] = { 0, 20 };

		double[][] demands = { { 7, 7, 7, 7, 7, 7 }, { 2, 3, 4, 5, 6, 7 }, { 8, 7, 6, 5, 4, 3 }, { 5, 6, 7, 8, 7, 6 },
				{ 8, 5, 2, 1, 2, 5 }, { 8, 4, 1, 3, 1, 3 }, { 1, 3, 8, 4, 8, 7 }, { 1, 4, 7, 3, 5, 8 },
				{ 3, 8, 4, 4, 6, 2 }, { 3, 1, 5, 8, 4, 4 } };

		for (int iK = 0; iK < K.length; iK++)
			for (int iv = 0; iv < v.length; iv++)
				for (int idemand = 0; idemand < demands.length; idemand++)
					for (int iP = 1; iP < P.length; iP++)
						for (int iCash = 0; iCash < initialCash.length; iCash++)
							for (int ib = 0; ib < b.length; ib++) {

								double[] meanDemand = demands[idemand];
								double fixOrderCost = K[iK];
								double variCost = v[iv];
								double holdingCost = 1;
								double price = P[iP];
								double iniCash = initialCash[iCash];
								double interestRate = b[ib];

								double truncationQuantile = 0.999;
								int stepSize = 1;
								double minInventoryState = 0;
								double maxInventoryState = 150;
								double minCashState = -100;
								double maxCashState = 800;
								double maxOrderQuantity = 50;
								double minCashRequired = -1000;

								// get demand possibilities for each period
								int T = meanDemand.length;
								Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
										.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other
																						// distributions
										.toArray(PoissonDist[]::new);
								double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

								// feasible actions
								Function<CashState, double[]> getFeasibleAction = s -> {
									double maxQ = (int) Math.min(maxOrderQuantity,
											Math.max(0, (s.getIniCash() - minCashRequired - fixOrderCost) / variCost));
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
									double cashBalanceBefore = state.getIniCash() + revenue - fixedCost - variableCost
											- holdCosts;
									double interest = interestRate * Math.max(-cashBalanceBefore, 0);
									double cashBalanceAfter = cashBalanceBefore - interest;
									double cashIncrement = cashBalanceAfter - state.getIniCash();
									return cashIncrement;
								};

								// state transition function
								StateTransitionFunction<CashState, Double, Double, CashState> stateTransition = (state,
										action, randomDemand) -> {
									double nextInventory = Math.max(0, state.getIniInventory() + action - randomDemand);
									double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
									double fixedCost = action > 0 ? fixOrderCost : 0;
									double variableCost = variCost * action;
									double holdCosts = holdingCost * Math.max(nextInventory, 0);
									double nextCash = state.getIniCash() + revenue - fixedCost - variableCost
											- holdCosts;
									nextCash = nextCash - Math.max(-nextCash, 0) * interestRate;
									nextCash = nextCash > maxCashState ? maxCashState : nextCash;
									nextCash = nextCash < minCashState ? minCashState : nextCash;
									nextInventory = nextInventory > maxInventoryState ? maxInventoryState
											: nextInventory;
									nextInventory = nextInventory < minInventoryState ? minInventoryState
											: nextInventory;
									// cash is integer or not
									nextCash = Math.round(nextCash * 10) / 10.0;
									return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
								};

								/*******************************************************************
								 * Solve
								 */
								CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction,
										stateTransition, immediateValue);
								int period = 1;
								double iniInventory = 0;
								CashState initialState = new CashState(period, iniInventory, iniCash);
								long currTime = System.currentTimeMillis();
								recursion.setTreeMapCacheAction();
								double finalValue = iniCash + recursion.getExpectedValue(initialState);
								System.out.println("final optimal cash is: " + finalValue);
								System.out.println("optimal order quantity in the first priod is : "
										+ recursion.getAction(initialState));
								double time = (System.currentTimeMillis() - currTime) / 1000;
								System.out.println("running time is " + time + "s");

								/*******************************************************************
								 * Find (s, S) and simulate
								 */
								int sampleNum = 10000;
								CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion);
								System.out.println("");
								double[][] optTable = recursion.getOptTable();
								FindsSOverDraft findsS = new FindsSOverDraft(T, iniCash);
								double[][] optsS = findsS.getsS(optTable);
								double simsSFinalValue = simuation.simulatesSOD(initialState, optsS, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
								double gap = (finalValue -simsSFinalValue)/finalValue*100;
								System.out.printf("Optimality gap is: %.2f%%\n", gap);
								
								String out = fixOrderCost + ",\t"+
										variCost + ",\t"+
										holdingCost + ",\t"+
										iniInventory + ",\t"+
										price + ",\t" +
										iniCash + ",\t" +
										interestRate + ",\t" +
										(idemand + 1) + ",\t" +
										finalValue + ",\t" +
										time + ",\t" + 
										simsSFinalValue + ",\t" +
										gap;	
								WriteToCsv.writeToFile("./" + "test_results.csv", out);
							}							
	}

}
