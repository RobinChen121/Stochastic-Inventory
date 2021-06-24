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
* @version 2018, April 22th, 9:38:22 
* @value testing for overdraft limit problem
*/
public class CashOverdraftLimitTesting{
	
	public static void main(String[] args) {
		String headString = "K, v, h, I0, pai, B0, minCash, rate, DemandPatt, OpValue, Time(sec), simsSValue, gap";
		WriteToCsv.writeToFile("./" + "test_results.csv", headString);

		double[] K = {10,15};
		double[] v = {1,2};
		double[] S = {5,10}; // selling price
		double salvageValue = 0.5;

		double[] minCash = {-40,-80};
		double[] b = {0.1,0.2};

		int initialCash[] = {0,20};
		double discountFactor = 0.95;

		double[][] demands = {{7,7,7,7,7,7},
				{2,3,4,5,6,7},
				{8,7,6,5,4,3},
				{5,6,7,8,7,6},
				{8,5,2,1,2,5},
				{8,4,1,3,1,3},
				{1,3,8,4,8,7},
				{1,4,7,3,5,8},
				{3,8,4,4,6,2},
				{3,1,5,8,4,4}
		};

		for (int iB = 0; iB < initialCash.length; iB++) 
			for (int iv = 0; iv < v.length; iv++) 
				for (int iK = 0; iK < K.length; iK++) 
					for (int is = 0; is < S.length; is++)
						for (int idemand = 0; idemand < demands.length; idemand++) 
							for (int ib = 0; ib < b.length; ib++)
								for (int iC = 0; iC < minCash.length; iC++){ 
		
		 
									double[] meanDemand = demands[idemand];
									double fixOrderCost = K[iK];
									double variCost = v[iv];
									double holdingCost = 1;
									double price =S[is];
									double iniCash = initialCash[iB];
									double interestRate = b[ib];
									double minCashRequired = minCash[iC];
									
									double truncationQuantile = 0.999;
									int stepSize = 1;
									double minInventoryState = 0;
									double maxInventoryState = 150;
									double minCashState = -100;
									double maxCashState = 800;
									double maxOrderQuantity = 50;							

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
											stateTransition, immediateValue, discountFactor);
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
									 * Find (s, C, S) and simulate
									 */
									int sampleNum = 10000;
									CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor);
									System.out.println("");
									double[][] optTable = recursion.getOptTable();
									FindsSOverDraft findsCS = new FindsSOverDraft(T, iniCash);
									double[][] optsCS = findsCS.getsCS(optTable);
									double simsSFinalValue = simuation.simulatesCSDraft(initialState, optsCS, minCashRequired, maxOrderQuantity, fixOrderCost, variCost);
									double gap = (finalValue -simsSFinalValue)/finalValue*100;
									System.out.printf("Optimality gap is: %.2f%%\n", gap);
									
									String out = fixOrderCost + ",\t"+
											variCost + ",\t"+
											holdingCost + ",\t"+
											iniInventory + ",\t"+
											price + ",\t" +
											minCashRequired + ",\t" +
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
