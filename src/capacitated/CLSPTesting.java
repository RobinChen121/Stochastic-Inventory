package capacitated;

import java.util.function.Function;
import java.util.stream.IntStream;

import sdp.inventory.GetPmf;
import sdp.inventory.Recursion;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.Recursion.OptDirection;
import sdp.inventory.Simulation;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToCsv;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
* @author Zhen Chen
* @date: 2018年11月14日 下午7:06:59  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  this is class to test clsp on large test bed. 
*/

public class CLSPTesting {

	public static void main(String[] args) {
		String headString = "K, v, h, I0, pai, coeVar, DemandPatt, OpValue, Time(sec), simValue";
		WriteToCsv.writeToFile("./" + "test_results.csv", headString);
		
		double[][] demands =
			 {{10,	10,	10,	10,	10,	10,	10,	10},
			 {15,	16,	15,	14,	11,	7,	6,	3},
			 {3,	6,	7,	11,	14,	15,	16,	15},
			 {15,   4,	4,	10,	18,	 4,	 4,	10},
			 {12,	7,	7,	10,	13,	7,	7,	12},
			 {2,	4,	7,	3,	10,	10,	3,	3},
			 {5,	15,	26,	44,	24,	15,	22,	10},
			 {4,	23,	28,	50,	39,	26,	19,	32},
			 {11,	14,	7,	11,	16,	31,	11,	48},
			 {18,	6,	22,	22,	51,	54,	22,	21},
		};
		
		double[] K = { 200, 300, 400 }; 
		double[] v = { 0, 1 };
		double[] pai = { 5, 10, 20};		
		double[] coeVar = {0.1, 0.2, 0.3};
		double holdingCost = 1;
		
		int maxOrderQuantity = 500;
		double minInventory = -500;
		double maxInventory = 500;
		double truncationQuantile = 0.9999;
		double stepSize = 1;
		
		for (int idemand = 0; idemand < demands.length; idemand++)
			for (int iv = 0; iv < v.length; iv++) 
				for (int ipai = 0; ipai < pai.length; ipai++) 
					for (int iK = 0; iK < K.length; iK++) 
						for (int iCoe = 0; iCoe < coeVar.length; iCoe++) {
							double[] meanDemand = demands[idemand];
							double fixedOrderingCost = K[iK];
							double variOrderingCost = v[iv];
							double penaltyCost = pai[ipai];
							double coeVarValue = coeVar[iCoe];
							
							int T = meanDemand.length;
							Distribution[] distributions = IntStream.iterate(0, i -> i + 1)
                                    .limit(T)
                                    //.mapToObj(i -> new PoissonDist(meanDemand[i]))
                                    .mapToObj(i -> new NormalDist(meanDemand[i], coeVarValue * meanDemand[i]))
                                    .toArray(Distribution[]::new);
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
								variableCost = variOrderingCost * action;
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
							int period = 1;
							double iniInventory = 0;
							State initialState = new State(period, iniInventory);
							double finalValue = recursion.getExpectedValue(initialState);
							long currTime = System.currentTimeMillis();
							System.out.println("final optimal expected value is: " + recursion.getExpectedValue(initialState));
							System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
							double time = (System.currentTimeMillis() - currTime) / 1000;
							System.out.println("running time is " + time + "s");
							
							/*******************************************************************
							 * Simulate
							 */
							int sampleNum = 10000;
							Simulation simuation = new Simulation(distributions, sampleNum, recursion);
							double simFinalValue = simuation.simulateSDPGivenSamplNum(initialState);
							System.out.println("***************************************************");
							
							String out = fixedOrderingCost + ",\t" 
									+ variOrderingCost + ",\t" 
									+ holdingCost + ",\t" 
									+ iniInventory + ",\t" 
									+ penaltyCost + ",\t" 
									+ coeVar[iCoe] + ",\t"
									+ (idemand + 1) + ",\t" 
									+ finalValue + ",\t" 
									+ time + ",\t" 
									+ simFinalValue;
							WriteToCsv.writeToFile("./" + "test_results.csv", out);
						}

	}

}


