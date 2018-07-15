package capacitated.fitss;

import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.IntStream;

import sdp.inventory.GetPmf;
import sdp.inventory.Recursion;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.Recursion.OptDirection;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToCsv;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 13, 2018---11:50:12 AM
*@description:  two level fitted (s, S) policy 
*/

public class TwoLevelFitsSTest {

	public static void main(String[] args) {
		String headString = "K, v, h, I0, pai, Qmax, DemandPatt, OpValue, Time(sec), simValue, error";
		WriteToCsv.writeToFile("./" + "test_results.csv", headString);

		double[][] demands = {{30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
				{6.6,9.3,11.1,12.9,16.8,21.6,24,26.4,31.5,33.9,36.3,40.8,44.4,47.1,48.3,50.1,50.1,48.9,47.1,44.4},
				{45.9,48.6,49.8,49.8,48.6,45.9,42.3,37.8,35.4,33,30.3,27.9,25.5,23.1,20.7,18.3,14.4,10.8,8.1,6.3},
				{36.3,30,23.7,21,23.7,30,36.3,39,36.3,30,23.7,21,23.7,30,36.3,30.9,24.3,21.3,26.4,33},
				{47.1,30,12.9,6,12.9,30,47.1,54,47.1,30,12.9,6,12.9,30,47.1,29.7,15,7.5,11.4,29.7},
				{62.7,27.3,9.9,23.7,1,22.8,32.7,34.5,67.2,6.6,14.4,40.8,3.9,63.3,25.5,45,53.1,24.6,10.2,49.8},
				{1,15.3,45.6,140.1,80.4,146.7,133.8,74.4,84.3,108.9,46.5,87.9,66,27.9,32.1,88.5,161.7,36,32.4,49.8},
				{14.1,24.3,70.8,118.2,49.2,86.1,152.4,117.3,226.2,208.2,78.3,58.5,96,33.3,57.3,115.5,17.7,134.7,127.8,179.7},
				{13.2,34.8,79.2,43.2,43.8,59.4,22.2,54.9,61.2,34.2,49.5,95.4,35.7,144.6,160.2,103.8,150.6,86.4,122.7,63.6},
				{14.7,56.4,19.2,83.7,135.9,67.2,66.9,155.1,87.3,164.1,193.8,67.2,64.5,132,34.8,131.4,132.9,35.7,172.5,152.1},
		};

		double[] K = {500,800,200};// 1500, 1000, 500
		double[] v = {0,5,10};
		double[] pai = {15, 10, 5}; // 15, 10, 5
		double[] capacity = {2, 3, 4}; // 3, 5, 7

		double truncationQuantile = 0.9999;  
		double stepSize = 1; 
		double minInventory = -300;
		double maxInventory = 800;
		double holdingCost = 1;

		for (int iK = 0; iK < K.length; iK++) {
			for (int iv = 0; iv < v.length; iv++) {
				for (int ipai = 0; ipai < pai.length; ipai++) {
					for (int idemand = 0; idemand < demands.length; idemand++) 
						for ( int icapacity = 0; icapacity < capacity.length; icapacity++){	      
							
							//double[] meanDemand = demands[idemand];
							double[] meanDemand = {9, 23, 53, 29};
							double fixedOrderingCost = K[iK];
							double proportionalOrderingCost = v[iv];
							double penaltyCost = pai[ipai];
							int maxOrderQuantity = (int) (Math
									.round(Arrays.stream(meanDemand).sum() / meanDemand.length) * capacity[icapacity]);
							
							// get demand possibilities for each period
							int T = meanDemand.length;
							Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
									.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other distributions
									.toArray(PoissonDist[]::new);
							double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

							// feasible actions
							Function<State, double[]> getFeasibleAction = s -> {
								double[] feasibleActions = new double[(int) (maxOrderQuantity / stepSize) + 1];
								int index = 0;
								for (double i = 0; i <= maxOrderQuantity; i = i + stepSize) {
									feasibleActions[index] = i;
									index++;
								}
								return feasibleActions;
							};

							// state transition function
							StateTransitionFunction<State, Double, Double, State> stateTransition = (state, action, randomDemand) -> {
								double nextInventory = state.getIniInventory() + action - randomDemand;
								nextInventory = nextInventory > maxInventory ? maxInventory : nextInventory;
								nextInventory = nextInventory < minInventory ? minInventory : nextInventory;
								return new State(state.getPeriod() + 1, nextInventory);
							};

							// immediate value
							ImmediateValueFunction<State, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
								double fixedCost = 0, variableCost = 0, inventoryLevel = 0, holdingCosts = 0, penaltyCosts = 0;
								fixedCost = action > 0 ? fixedOrderingCost : 0;
								variableCost = proportionalOrderingCost * action;
								inventoryLevel = state.getIniInventory() + action - randomDemand;
								holdingCosts = holdingCost * Math.max(inventoryLevel, 0);
								penaltyCosts = penaltyCost * Math.max(-inventoryLevel, 0);
								double totalCosts = fixedCost + variableCost + holdingCosts + penaltyCosts;
								return totalCosts;
							};

							/*******************************************************************
							 * Solve
							 */
							Recursion recursion = new Recursion(OptDirection.MIN, pmf, getFeasibleAction, stateTransition, immediateValue);
							recursion.setTreeMapCacheAction(); // using tree map when finding levels for s and S
							int period = 1;
							double iniInventory = 0;
							State initialState = new State(period, iniInventory);
							long currTime = System.currentTimeMillis();
							double finalValue = recursion.getExpectedValue(initialState);
							System.out.println("final optimal expected value is: " + finalValue);
							System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
							double time = (System.currentTimeMillis() - currTime) / 1000;
							System.out.println("running time is " + time + "s");

							/*******************************************************************
							 * Simulating two level s S results
							 */
							int sampleNum = 10000;
							SimulateFitsS simuation = new SimulateFitsS(distributions, sampleNum, recursion);
							double[][] optTable = recursion.getOptTable();
							FindsS findsS = new FindsS(maxOrderQuantity, T);
							double[][] optsS = findsS.getTwosS(optTable);
							System.out.println(Arrays.deepToString(optsS));
							simuation.simulateSDPGivenSamplNum(initialState);
							double simFinalValue = simuation.simulateTwosS(initialState, optsS, maxOrderQuantity);
							System.out.printf("Optimality gap is: %.2f%%\n",
									(simFinalValue - finalValue) / finalValue * 100);
							String out = fixedOrderingCost + ",\t" 
									+ proportionalOrderingCost + ",\t" 
									+ holdingCost + ",\t" 
									+ iniInventory + ",\t" 
									+ penaltyCost + ",\t" 
									+ maxOrderQuantity + ",\t"
									+ (idemand + 1) + ",\t" 
									+ finalValue + ",\t" 
									+ time + ",\t" 
									+ simFinalValue + ",\t"
									+ (simFinalValue - finalValue) / finalValue;
							WriteToCsv.writeToFile("./" + "test_results.csv", out);
						}
				}
			}
		}
	}
}

