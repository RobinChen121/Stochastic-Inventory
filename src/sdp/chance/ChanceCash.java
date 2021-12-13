package sdp.chance;


import java.awt.print.Printable;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.function.IntPredicate;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.sun.net.httpserver.Authenticator.Result;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import milp.LostSaleChance;
import milp.PositiveCashChance;
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
import umontreal.ssj.probdist.LognormalDist;
import umontreal.ssj.probdist.PoissonDist;


/**
 * @author chen
 * @email: okchen321@163.com
 * @date: 2021 Jun 18, 17:41:22  
 * @desp: try the chance programming approach for the stochastic cash-flow inventory problem.
 * 
 * 
 * objective is to maximize the survival probability.
 * 
 * big scenario number will increase the performance of extended SAA.
 * 
 * sample numbers [5, 5, 3, 3, 3, 3] are 2025s for saa, extended saa 3747s, service rate is 60% 
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
		double iniCash = 150;
		double iniI = 0;
		double trunQuantile = 0.9999;
		double serviceRate = 0.6; // maximum negative possibility rate is 1 - serviceRate. the higher value results in slower running speed.

		int[] sampleNums = {5, 5, 5, 3, 3, 3}; // sample number in each period, the number of samples in the first period can have big influence
		double[] meanDemand = {10, 5, 5, 5, 5, 5};
		int T = sampleNums.length;
		
		double holdCostUnit = 0;
		double salvageValueUnit = 5;
		
		double[] prices = new double[T];
		double[] variCostUnits = new double[T];
		double[] overheadCosts = new double[T];		
//		double[] mus = new double[T];
//		double[] sigmas = new double[T];
		
		Arrays.fill(prices, 23);	
		Arrays.fill(variCostUnits, 10);
		Arrays.fill(overheadCosts, 100); // overhead costs
//		Arrays.fill(mus, 3.6);
//		Arrays.fill(sigmas, 0.6);
		
		int sampleNum = 200;  // simulating sample number in testing SAA and extended SAA
		
		double maxOrderQuantity = 300; // maximum ordering quantity when having enough cash

		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				// .mapToObj(i -> new NormalDist(meanDemand[i], Math.sqrt(meanDemand[i]))) //can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i])).toArray(Distribution[]::new); // Poisson demand
				//.mapToObj(i -> new LognormalDist(mus[i], sigmas[i])).toArray(Distribution[]::new);

		/**
		 * solve the problem by SAA
		 */
		// generate scenarios, samples in each period form a scenario tree
		Sampling sampling = new Sampling();
		double[][] scenarios = sampling.generateLHSamples(distributions, sampleNums);
		
		int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);		
		int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
		
		LostSaleChance model = new LostSaleChance(distributions, sampleNums, iniCash, iniI, prices, variCostUnits, 
							salvageValueUnit, holdCostUnit, overheadCosts, serviceRate, scenarios);
		long currTime = System.currentTimeMillis();
		double[] result;
		double time1;
		double positiveScenario;
		double survivalProb;
		double lostRate;
		NumberFormat nf = NumberFormat.getPercentInstance();
		nf.setMinimumFractionDigits(5);		
		DecimalFormat df = new DecimalFormat("###, ###");
		
		result = model.solveMaxSurvival();		
		
		time1 = (System.currentTimeMillis() - currTime) / 1000.00;
	    currTime = System.currentTimeMillis();	    
	    System.out.println("**********************************************");
	    System.out.println("result of SAA-scenario tree before sorting scenarios: ");
	    System.out.println("running time is " + time1 + "s");	
	    System.out.printf("first stage decison Q is: %.2f\n", result[0]);
	    positiveScenario = result[1];
	    System.out.printf("Objective value is: %.0f in %d scenarios\n", result[1], sampleNumTotal);
	    survivalProb = 100 * result[1] / sampleNumTotal;
	    System.out.printf("Survival probability is: %.5f%%\n", survivalProb);
	    System.out.println("lost sale scenario number in the solution is : " + result[2]);
	    System.out.println("maximum lost sale scenario number allowed is: " + negativeScenarioNumRequire);
	    lostRate = result[2] / (double) sampleNumTotal;
	    System.out.println("lost sale rate of SAA is: " + nf.format(lostRate));
	    System.out.println("lost sale max required rate is: " + nf.format(1 - serviceRate));
	    System.out.println();
	    
	    /**
		 * Simulate the restult of SAA
		 */
	    int stepSize = 1;
	    double fixOrderCost = 0;
	    double depositeRate = 0;
	    double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -1000; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;
		double discountFactor = 1;
		double[][][] pmf = new GetPmf(distributions, trunQuantile, stepSize).getpmf();	
		
		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			int t = state.getPeriod() - 1;
			double revenue = prices[t] * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCostUnits[t] * action;
			double deposite = (state.getIniCash() - fixedCost - variableCost) * (1 + depositeRate);
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdCostUnit * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue + deposite - holdCosts - overheadCosts[t] - state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValueUnit * Math.max(inventoryLevel, 0) : 0;
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
			nextCash = Math.round(nextCash * 1) / 1; // the right should be a decimal
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};
		
	    
		int period = 1;		
	    CashState initialState = new CashState(period, iniI, iniCash);
	    CashSimulation simulation1 = new CashSimulation(distributions, sampleNum, immediateValue, stateTransition); // no need to add overheadCost in this class
	    // result[0] = 14;
//	    double[] result1 = simulation1.simulateSAALostRate(initialState, result[0], serviceRate, sampleNums, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenarios, sampleNum);
//		System.out.println("final simulated survival probability of SAA in " + df.format(sampleNum) + " samples is: " + nf.format(result1[0]));
//		System.out.println("final simulated lost sale rate of SAA " + " is: " + nf.format(result1[1])); 
	    
		/**
		 * solve the problem by extended formulation of SAA
		 * 
		 */
		currTime = System.currentTimeMillis();
		result = model.solveSort();	
		
		time1 = (System.currentTimeMillis() - currTime) / 1000.00;
	    currTime = System.currentTimeMillis();	    
	    System.out.println("**********************************************");
	    System.out.println("after sorting, result of SAA-scenario tree: ");
	    System.out.println("running time is " + time1 + "s");	
	    System.out.printf("first stage decison Q is: %.2f\n", result[0]);
	    positiveScenario = result[1];
	    System.out.printf("Objective value is: %.0f in %d scenarios\n", positiveScenario, sampleNumTotal);
	    survivalProb = 100 * result[1] / sampleNumTotal;
	    System.out.printf("Survival probability is: %.5f%%\n", survivalProb);
	    System.out.println("lost sale scenario number in the solution is : " + result[2]);
	    System.out.println("maximum lost sale scenario number allowed is: " + negativeScenarioNumRequire);
	    lostRate = result[2] / (double) sampleNumTotal;
	    System.out.println("lost sale rate of SAA is: " + nf.format(lostRate));
	    System.out.println("lost sale max required rate is: " + nf.format(1 - serviceRate));
	    System.out.println();
	    
	    /**
		 * Simulate the restult of extended SAA
		 */	    
//	    result1 = simulation1.simulateExtendSAALostRate(initialState, result[0], serviceRate, sampleNums, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenarios, sampleNum);
//		System.out.println("final simulated survival probability of extended SAA in " + df.format(sampleNum) + " samples is: " + nf.format(result1[0]));
//		System.out.println("final simulated lost sale rate of extended SAA " + " is: " + nf.format(result1[1])); 	    
		
		
	    /**
		 * solve the problem by SDP to get a lower bound
		 */	
		// feasible actions
		// in fact, no cash constraint in this paper
		Function<CashState, double[]> getFeasibleAction1 = s -> {
			int t = s.getPeriod() - 1;
			double thisPeriodServRate = (1 - serviceRate) / T;
//			double minQ = Math.ceil(distributions[t].inverseF(1 - thisPeriodServRate)); // minimum ordering quantity in each period 
			double minQ = 0;
			return DoubleStream.iterate(minQ, i -> i + stepSize).limit((int) maxOrderQuantity + 1).toArray();
		};
		
		
		/*******************************************************************
		 * Solve
		 */
		CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction1, stateTransition,
				immediateValue, discountFactor);
		currTime = System.currentTimeMillis();
		recursion.setTreeMapCacheAction();
		double finalValue = recursion.getSurvProb(initialState);
		System.out.println("**********************************************");
		System.out.println("result of SDP with no service rate constraint is: ");
		System.out.println("survival probability for this initial state is: " + nf.format(finalValue));
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		/*******************************************************************
		 * Simulating sdp results
		 */		
		sampleNum = 10000;
		CashSimulation simulation = new CashSimulation(distributions, sampleNum, recursion, discountFactor); // no need to add overheadCost in this class
		double[] result2 = simulation.simulateSDPGivenSamplNumLostRate(initialState, immediateValue);
		System.out.println("final simulated survival probability in " + df.format(sampleNum) + " samples is: " + nf.format(result2[0]));
		System.out.println("final simulated lost sale rate " + " is: " + nf.format(result2[1])); 
		
		// feasible actions 2
		// in fact, no cash constraint in this paper
		Function<CashState, double[]> getFeasibleAction2 = s -> {
			int t = s.getPeriod() - 1;
			double thisPeriodServRate = (1 - serviceRate) / T;
			double minQ = Math.ceil(distributions[t].inverseF(1 - thisPeriodServRate)); // minimum ordering quantity in each period 
			return DoubleStream.iterate(minQ, i -> i + stepSize).limit((int) maxOrderQuantity + 1).toArray();
		};
		
		/*******************************************************************
		 * Solve
		 */
		recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction2, stateTransition,
				immediateValue, discountFactor);
		period = 1;		
		initialState = new CashState(period, iniI, iniCash);
		currTime = System.currentTimeMillis();
		recursion.setTreeMapCacheAction();
		finalValue = recursion.getSurvProb(initialState);
		System.out.println("**********************************************");
		System.out.println("result of SDP with service rate constraint is: ");
		finalValue = finalValue * 100;
		System.out.printf("survival probability for this initial state is: %.2f%%\n", finalValue);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		/*******************************************************************
		 * Simulating sdp results
		 */		
		sampleNum = 10000;
		simulation = new CashSimulation(distributions, sampleNum, recursion, discountFactor); // no need to add overheadCost in this class
		result2 = simulation.simulateSDPGivenSamplNumLostRate(initialState, immediateValue);
		System.out.println("final simulated survival probability in " + df.format(sampleNum) + " samples is: " + nf.format(result2[0]));
		System.out.println("final simulated lost sale rate " + " is: " + nf.format(result2[1]));	
		
		
		
		
		
		
	}
}
