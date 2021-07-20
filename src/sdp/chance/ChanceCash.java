package sdp.chance;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.function.IntPredicate;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import milp.GurobiChance;
import milp.SimulateChanceCash;
import sdp.cash.CashRecursion;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.CartesianProduct;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;


/**
 * @author chen
 * @email: okchen321@163.com
 * @date: 2021 Jun 18, 17:41:22  
 * @desp: try the chance programming approach for the stochastic inventory problem.
 *
 */
public class ChanceCash {
	
	static int[][] scenarioIndexes(int[] sampleNums, int sampleNumTotal){
		int T = sampleNums.length;
		int[][] arr = new int[sampleNumTotal][T];
		for (int i = 0; i < sampleNumTotal; i++) {
			for (int t = 0; t < T; t++) {
				arr[i][t] = 0;
			}
		}
		return null;
	}
	
	public static void main(String[] args) {
		double iniCash = 10;
		double iniI = 0;
		double variCostUnit = 1;
		double salvageValueUnit = 0.5;
		double trunQuantile = 0.9999;
		double serviceRate = 0; // maximum negative possibility rate is 1 - serviceRate

		double[] price = {2, 2, 2, 2};
		double[] meanDemand = {50, 50, 50};
		int[] sampleNums = {10, 10, 10}; // sample number in each period
		double maxOrderQuantity = 150; // maximum ordering quantity when having enough cash
		int T = sampleNums.length;
		double[] overheadCost = new double[T];
		Arrays.fill(overheadCost, 5); // overhead costs
		

		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				// .mapToObj(i -> new NormalDist(meanDemand[i], Math.sqrt(meanDemand[i]))) //can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i])).toArray(Distribution[]::new); // Poisson demand

		/**
		 * solve the problem by SAA
		 */
		// generate scenarios, samples in each period form a scenario tree
		Sampling sampling = new Sampling();
		double[][] scenarios = sampling.generateLHSamples(distributions, sampleNums);
		int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);		
		int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
		
	    System.out.println("result of SAA-scenario tree: ");
		GurobiChance model = new GurobiChance(distributions, sampleNums, iniCash, iniI, price, variCostUnit, 
							salvageValueUnit, overheadCost, serviceRate, scenarios);
		long currTime = System.currentTimeMillis();
		double[] result = model.solve();		
		double time1 = (System.currentTimeMillis() - currTime) / 1000.00;
	    currTime = System.currentTimeMillis();
	    double[] result2 = model.solveSort();	    
	    System.out.println();
	    System.out.println("**********************************************");
	    System.out.println("result of SAA-scenario tree before sorting scenarios: ");
	    System.out.println("running time is " + time1 + "s");	
	    System.out.printf("first stage decison Q is: %.2f\n", result[0]);
	    System.out.printf("Objective value is: %.2f\n", result[1]);
	    System.out.println("negative scenario number in the solution is : " + result[2]);
	    System.out.println("maximum negative scenario number required is: " + negativeScenarioNumRequire);
	    
	    System.out.println();
		System.out.println("**********************************************");
	    System.out.println("result of SAA-scenario tree after sorting scenarios: ");
		currTime = System.currentTimeMillis();	
		double time3 = (System.currentTimeMillis() - currTime) / 1000.00;
	    currTime = System.currentTimeMillis();  

	    System.out.println("running time is " + time3 + "s");	
	    System.out.printf("first stage decison Q is: %.2f\n", result2[0]);
	    System.out.printf("Objective value is: %.2f\n", result2[1]);
	    System.out.println("negative scenario number in the solution is : " + result2[2]);
	    System.out.println("maximum negative scenario number required is: " + negativeScenarioNumRequire);
	    
//	    System.out.println();
//	    System.out.println("**********************************************");
//	    System.out.println("result of SAA-No-scenario tree: ");
//		double time2 = (System.currentTimeMillis() - currTime) / 1000.00;
//		System.out.println("running time is " + time2 + "s");	
//	    System.out.printf("first stage decison Q is: %.2f\n", result2[0]);
//	    System.out.printf("Objective value is: %.2f\n", result2[1]);
	    
	    /**
		 * solve the problem by SDP to get a lower bound
		 */
	    int stepSize = 1;
	    double fixOrderCost = 0;
	    double depositeRate = 0;
	    double holdingCost = 0;
	    double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -1000; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;
		double discountFactor = 1;
		double[][][] pmf = new GetPmf(distributions, trunQuantile, stepSize).getpmf();		
		
		// feasible actions
		Function<CashState, double[]> getFeasibleAction = s -> {
			int t = s.getPeriod() - 1;
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (s.getIniCash() - fixOrderCost) / variCostUnit));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			int t = state.getPeriod() - 1;
			double revenue = price[t] * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCostUnit * action;
			double deposite = (state.getIniCash() - fixedCost - variableCost) * (1 + depositeRate);
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue + deposite - holdCosts - overheadCost[t] - state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValueUnit * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			//double endCash = state.getIniCash() + cashIncrement;
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
			nextCash = Math.round(nextCash * 10) / 10.0; // the right should be a decimal
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};

		/*******************************************************************
		 * Solve
		 */
		CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition,
				immediateValue, discountFactor);
		int period = 1;		
		CashState initialState = new CashState(period, iniI, iniCash);
		currTime = System.currentTimeMillis();
		recursion.setTreeMapCacheAction();
		double finalValue = iniCash + recursion.getExpectedValue(initialState);
		System.out.println();
	    System.out.println("**********************************************");
	    System.out.println("result of SDP: ");
	    double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		System.out.println("first stage decision Q is : " + recursion.getAction(initialState));
		System.out.println("final optimal cash  is " + finalValue);    
		
		/*******************************************************************
		 * Simulating sdp results
		 */
		int sampleNum = 10000;		
		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor); // no need to add overheadCost in this class
		double simFinalValue = simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001; 
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);	    
		
		/*******************************************************************
		 * Simulating SAA result
		 */
		SimulateChanceCash simulate = new SimulateChanceCash(distributions, iniCash, iniI, price, variCostUnit, salvageValueUnit, overheadCost, serviceRate, stateTransition, 
				immediateValue, discountFactor, sampleNums, scenarios);
		double simFinalValue2 = simulate.simulateSDPGivenSamplNum(initialState, result[0], 1000);
		System.out.println("**********************************************");
		System.out.println("simulate result for chanced SAA is " + simFinalValue2);  
	}
}
