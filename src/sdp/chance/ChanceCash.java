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
import sdp.cash.RiskSimulation;
import sdp.cash.CashState;
import sdp.cash.RiskState;
import sdp.cash.RiskRecursion;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.CartesianProduct;
import sdp.sampling.Sampling;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.LognormalDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;


/**
 * @author chen
 * @email: okchen321@163.com
 * @date: 2021 Jun 18, 17:41:22  
 * @desp: try the chance programming approach for the stochastic cash-flow inventory problem.
 * 
 * 
 * objective is to maximize the survival probability.
 * large scenario number will improve the performance of extended SAA.
 * when there is inventory holding cost, sometimes the results are strange: lost sale rate in the simulation too large.
 *
 */
public class ChanceCash {
	
	static int[][] scenarioIndexes(int[] sampleNumSims, int sampleNumSimTotal){
		int T = sampleNumSims.length;
		int[][] arr = new int[sampleNumSimTotal][T];
		for (int i = 0; i < sampleNumSimTotal; i++) {
			for (int t = 0; t < T; t++) {
				arr[i][t] = 0;
			}
		}
		return null;
	}
	
	public static void main(String[] args) {
		
		WriteToExcelTxt wr = new WriteToExcelTxt();
		String fileName = "nita_N.xls";
		String headString =  
				"sampleNumSim" + "\t" +  "etaFrac" + "\t" + "eta" + "\t" + "SAAObj" + "\t" + "simSAAObj" + "\t" +
		         "simSAAService" + "\t" + "scenarioObj" + "\t" + "simScenarioObj" + "\t" + "simScenarioService"
		         + "\t" + "SDPObj" + "\t" + "simSDPService" + "\t" + "SDPLbObj" + "\t" + "simSDPLbService" + "\t"
		         + "rollingObj" + "\t" + "rollingService" + "\t" + "rollingLength";
		wr.writeToFile(fileName, headString);
		
		
		double[] meanDemand = {11, 14, 18, 23, 28, 33, 38, 42, 46, 49, 50, 49};
		double iniCash = 130;
		int T = meanDemand.length;
		
		double minInventoryState = 0;
		double maxInventoryState = 300;
		double minCashState = -10; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 200;
		double discountFactor = 1;		
		
		double iniI = 0;
		double trunQuantile = 0.9999;
		
		double serviceRate = 0.95; // the higher value results in slower running speed. maximum negative possibility rate is 1 - serviceRate. 
		int etaFrac = 1;
		double serviceRateSAA = 1 - (1- serviceRate)/etaFrac;
			
		int sampleNumPeriod = 10;
		int[] sampleNumSims = new int[T]; // sample number in each period, the number of samples in the first period can have big influence		
		int[] sampleNumSimsRolling = new int[T];
		Arrays.fill(sampleNumSimsRolling, sampleNumPeriod);
		Arrays.fill(sampleNumSims, sampleNumPeriod);
		
		int rollingLength = 2; // rolling horizon length
		double meanDemandSum = Arrays.stream(meanDemand).sum();
		double rollingDemandSum = Arrays.stream(meanDemand).limit(rollingLength).sum();
		double portion = rollingDemandSum / meanDemandSum;
		double rollingServiceRate = Math.pow(serviceRateSAA, portion);		
		int sampleNumSim = 100;  // simulating sample number in testing SAA and extended SAA
		
		double holdCostUnit = 0;
		double salvageValueUnit = 0.5;	
		
		double[] prices = new double[T];
		double[] variCostUnits = new double[T];
		double[] overheadCosts = new double[T];		
		
		Arrays.fill(prices, 5);	
		Arrays.fill(variCostUnits, 1);
		Arrays.fill(overheadCosts, 80); // overhead costs
		
		
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash
	    int stepSize = 1;
		double sigmaCoe = 0.25;
		
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				//.mapToObj(i -> new NormalDist(meanDemand[i], sigmaCoe * meanDemand[i])).toArray(Distribution[]::new);// can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i])).toArray(Distribution[]::new);
		
//		for (int runTime = 0; runTime < 10; runTime++) {
		
		
		/**
		 * solve the problem by SAA
		 */
		// generate samples
		Sampling sampling = new Sampling();
		double[][] scenarios = sampling.generateLHSamples(distributions, sampleNumSims);
	
		int sampleNumSimTotal = IntStream.of(sampleNumSims).reduce(1, (a, b) -> a * b);		
		int sampleNumSimTotalSimulate = IntStream.of(sampleNumSimsRolling).limit(rollingLength).reduce(1, (a, b) -> a * b);
		int negativeScenarioNumRequire = (int) (sampleNumSimTotal * (1 - serviceRateSAA));
		

		long currTime = System.currentTimeMillis();
		double[] result;
		double time1;
		double positiveScenario;
		double ObjSAA;
		double lostRate;
		NumberFormat nf = NumberFormat.getPercentInstance();
		nf.setMinimumFractionDigits(5);		
		DecimalFormat df = new DecimalFormat("###, ###");
	    double fixOrderCost = 0;
	    double depositeRate = 0;
		double[][][] pmf = new GetPmf(distributions, trunQuantile, stepSize).getpmf();	
		int period = 1;		
	    RiskState initialState = new RiskState(period, iniI, iniCash, false);
	    double[] resultSim; 
	    
		// immediate value
		ImmediateValueFunction<RiskState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			int t = state.getPeriod() - 1;
			double revenue = prices[t] * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCostUnits[t] * action;
			double deposite = (state.getIniCash() - fixedCost - variableCost) * (1 + depositeRate);
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdCostUnit * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue + deposite - holdCosts - overheadCosts[t]
					- state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValueUnit * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			return cashIncrement;
		};
		
			// state transition function
		StateTransitionFunction<RiskState, Double, Double, RiskState> stateTransition = (state, action,
				randomDemand) -> {
			double nextInventory = Math.max(0, state.getIniInventory() + action - randomDemand);
			double nextCash = state.getIniCash() + immediateValue.apply(state, action, randomDemand);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			// cash is integer or not
			nextCash = Math.round(nextCash * 1) / 1; // the right should be a decimal
			boolean bankruptBefore = false;
			if (nextCash < 0)
				bankruptBefore = true;
			return new RiskState(state.getPeriod() + 1, nextInventory, nextCash, bankruptBefore);
		};
		
	    RiskSimulation simulation1 = new RiskSimulation(distributions, sampleNumSim, immediateValue, stateTransition); // no need to add overheadCost in this class
	    double error;


//		LostSaleChance model = new LostSaleChance(distributions, sampleNumSims, iniCash, iniI, prices, variCostUnits, 
//				salvageValueUnit, holdCostUnit, overheadCosts, serviceRateSAA, scenarios);
//        result = model.solveMaxSurvival();				
//		time1 = (System.currentTimeMillis() - currTime) / 1000.00;
//	    currTime = System.currentTimeMillis();	    
//	    System.out.println("*********************************************************************");
//	    System.out.println("result of SAA: ");
//	    System.out.println("running time is " + time1 + "s");	
//	    System.out.printf("first stage decison Q is: %.2f\n", result[0]);
//	    positiveScenario = result[1];
//	    System.out.printf("Objective value is: %.0f in %d scenarios\n", result[1], sampleNumSimTotal);
//	    ObjSAA = 100 * result[1] / sampleNumSimTotal;
//	    System.out.printf("Survival probability is: %.5f%%\n", ObjSAA);
//	    System.out.println("lost sale scenario number in the solution is : " + result[2]);
//	    System.out.println("maximum lost sale scenario number allowed is: " + negativeScenarioNumRequire);
//	    lostRate = result[2] / (double) sampleNumSimTotal;
//	    System.out.println("lost sale rate of SAA is: " + nf.format(lostRate));
//	    System.out.println("lost sale max required rate is: " + nf.format(1 - serviceRate));
//	    System.out.println();
	    
	    /**
		 * Simulate the restult of SAA
		 */		    
//	    double SAAServiceSim;	    
//	    currTime = System.currentTimeMillis();
//	    resultSim = simulation1.simulateSAA(initialState, result[0], serviceRateSAA, sampleNumSims, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenarios, sampleNumSim);
//	    time1 = (System.currentTimeMillis() - currTime) / 1000.00;  
//	    System.out.println("running time is " + time1 + "s");
//	    double objSAASim = resultSim[0];
//	    System.out.println("final simulated survival probability of SAA in " + df.format(sampleNumSim) + " samples is: " + nf.format(objSAASim));
//		error  = 1.96 * Math.sqrt(resultSim[1]*(1 - resultSim[1]) / sampleNumSim);
//		SAAServiceSim = 1-resultSim[1];
//		System.out.println("final simulated service sale rate of SAA " + " is: " + nf.format(SAAServiceSim) + " with error " + nf.format(error)); 
//	    
//		/**
//		 * solve the problem by scenario tree
//		 */
//		scenarios = sampling.generateLHSamples(distributions, sampleNumSims); // new scenarios
//		model = new LostSaleChance(distributions, sampleNumSims, iniCash, iniI, prices, variCostUnits, 
//				salvageValueUnit, holdCostUnit, overheadCosts, serviceRateSAA, scenarios);
//		currTime = System.currentTimeMillis();
//		result = model.solveScenario();	// same result with soveSort or solveSort2, but less computational time
//		                                   // former name is solveSortFurther()
		
//		time1 = (System.currentTimeMillis() - currTime) / 1000.00;    
//	    System.out.println("**************************************************************");
//	    System.out.println("result of scenario tree: ");
//	    System.out.println("running time is " + time1 + "s");	
//	    System.out.printf("first stage decison Q is: %.2f\n", result[0]);
//	    positiveScenario = result[1];
//	    System.out.printf("Objective value is: %.0f in %d scenarios\n", positiveScenario, sampleNumSimTotal);
//	    double ObjScenario = 100 * result[1] / sampleNumSimTotal;
//	    System.out.printf("Survival probability is: %.5f%%\n", ObjScenario);
//	    System.out.println("lost sale scenario number in the solution is : " + result[2]);
//	    System.out.println("maximum lost sale scenario number allowed is: " + negativeScenarioNumRequire);
//	    lostRate = result[2] / (double) sampleNumSimTotal;
//	    System.out.println("lost sale rate of scenario tree model is: " + nf.format(lostRate));
//	    System.out.println("lost sale max required rate is: " + nf.format(1 - serviceRate));
//	    System.out.println();
//	    
//	    /**
//		 * Simulate the result of scenario tree
//		 */	   
//	    currTime = System.currentTimeMillis();
//	    resultSim = simulation1.simulateScenarioTree(initialState, result[0], serviceRateSAA, sampleNumSims, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenarios, sampleNumSim);
//	    time1 = (System.currentTimeMillis() - currTime) / 1000.00;  
//	    double ObjScenarioSim = resultSim[0];
//	    System.out.println("running time is " + time1 + "s");
//	    System.out.println("final simulated survival probability of scenario tree model in " + df.format(sampleNumSim) + " samples is: " + nf.format(resultSim[0]));
//		error  = 1.96 * Math.sqrt(resultSim[1]*(1 - resultSim[1]) / sampleNumSim);
//		double ScenarioServiceSim = 1 - resultSim[1];
//		System.out.println("final simulated service rate of scenario tree model " + " is: " + nf.format(ScenarioServiceSim) + " with error " + nf.format(error));		
//
//	    /**
//		 * solve the problem by SDP when there is no joint chance constraint
//		 */	
//		// feasible actions
//		Function<RiskState, double[]> getFeasibleAction = s -> {
//			int t = s.getPeriod() - 1;
//			double maxQ = Math.min(s.iniCash/variCostUnits[t], maxOrderQuantity);
//			if (s.getBankruptBefore() == true)
//				maxQ = 0;
//			maxQ = Math.max(maxQ, 0);
//			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
//		};	
//		
//		/*******************************************************************
//		 * Solve
//		 */
//		RiskRecursion recursion = new RiskRecursion(pmf, getFeasibleAction, stateTransition, immediateValue);
//		double SDPSurvival;
//		double time;
//		
//		currTime = System.currentTimeMillis();
//		recursion.setTreeMapCacheAction();
//		SDPSurvival = recursion.getSurvProb(initialState);
//		System.out.println("**********************************************");
//		System.out.println("result of SDP is: ");
//		System.out.println("survival probability for this initial state is: " + nf.format(SDPSurvival));
//		double optQ = recursion.getAction(initialState);
//		System.out.println("optimal order quantity in the first priod is : " + optQ);
//		time = (System.currentTimeMillis() - currTime) / 1000;
//		System.out.println("running time is " + time + "s");
//		
//		/*******************************************************************
//		 * Simulating sdp results
//		 */	
//		sampleNumSim = 1000;
//		RiskSimulation simulation = new RiskSimulation(distributions, sampleNumSim, recursion); // no need to add overheadCost in this class
//		double[] result2 = simulation.simulateLostSale(initialState, immediateValue);
//		System.out.println("final simulated survival probability in " + df.format(sampleNumSim) + " samples is: " + nf.format(result2[0]));
//		System.out.println("final simulated lost sale rate " + " is: " + nf.format(result2[1])); 
//		double serviceSDP = 1 - result2[1];
//		System.out.println("final simulated service sale rate " + " is: " + nf.format(serviceSDP)); 
//		
//		/*******************************************************************
//		 * solve the problem by SDP when there is individual chance constraint approximation
//		 */			
//		// feasible actions 2
//		// in fact, no cash constraint in this paper
//		Function<RiskState, double[]> getFeasibleAction2 = s -> {
//			int t = s.getPeriod() - 1;
//			double thisPeriodServRate = 1- (1 - serviceRate) / T;
//			double minQ = Math.ceil(distributions[t].inverseF(thisPeriodServRate)); // minimum ordering quantity in each period 
//			double maxQ = Math.min(s.iniCash/variCostUnits[t], maxOrderQuantity);
//			if (s.getBankruptBefore() == true)
//				maxQ = 0;
//			if (maxQ < minQ) {
//				maxQ = 0;
//				minQ = 0;
//			}
//			maxQ = Math.max(maxQ, 0);
//			return DoubleStream.iterate(minQ, i -> i + stepSize).limit((int) maxQ + 1).toArray();
//		};
//		
//		/*******************************************************************
//		 * Solve
//		 */
//		recursion = new RiskRecursion(pmf, getFeasibleAction2, stateTransition, immediateValue);
//		period = 1;		
//		initialState = new RiskState(period, iniI, iniCash, false);
//		currTime = System.currentTimeMillis();
//		recursion.setTreeMapCacheAction();
//		double SDPLbSurvival;
//		SDPLbSurvival= recursion.getSurvProb(initialState);
//		System.out.println("**********************************************");
//		System.out.println("result of SDP with service rate constraint is: ");
//		System.out.println("survival probability for this initial state is: " + nf.format(SDPLbSurvival));
//		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
//		time = (System.currentTimeMillis() - currTime) / 1000;
//		System.out.println("running time is " + time + "s");
//		
//		/*******************************************************************
//		 * Simulating sdp results
//		 */		
//		sampleNumSim = 10000;
//		simulation = new RiskSimulation(distributions, sampleNumSim, recursion); // no need to add overheadCost in this class
//		result2 = simulation.simulateLostSale(initialState, immediateValue);
//		System.out.println("final simulated survival probability in " + df.format(sampleNumSim) + " samples is: " + nf.format(result2[0]));
//		System.out.println("final simulated lost sale rate " + " is: " + nf.format(result2[1]));
//		double serviceSDPLB = 1 - result2[1];
//		System.out.println("final simulated service rate " + " is: " + nf.format(serviceSDPLB));
		
		/**
		 * solve the problem by rolling horizon of SAA
		 * 
		 */
		sampleNumSim = 100; // number of scenarios for rolling SAA
	    currTime = System.currentTimeMillis();
	    System.out.println("**********************************************");
	    scenarios = sampling.generateLHSamples(distributions, sampleNumSimsRolling);
	    resultSim = simulation1.rollingHoirzon(rollingLength, initialState, rollingServiceRate, sampleNumSimsRolling, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenarios, sampleNumSim);
	    time1 = (System.currentTimeMillis() - currTime) / 1000.00;	    
	    System.out.println("after rolling horizon for length " + rollingLength +", " + "total horizon length is " + T + ", result is: ");
	    System.out.println("running time is " + time1 + "s");
	    double rollingObj = resultSim[0];
	    System.out.println("final simulated survival probability of rolling SAA in " + sampleNumSim + " samples is: " + nf.format(rollingObj));
	    double sigma2 = Math.sqrt(resultSim[1]*(1 - resultSim[1])/sampleNumSim);
		double error2  = 1.96*sigma2;
		double serviceRateRolling = 1 - resultSim[1];
		System.out.println("the service rate for simulated SAA rolling horizon is " + nf.format(serviceRateRolling) + ", with error " + nf.format(error2));
						
		/**
		 * output to excel
		 * 
		 */
//		double[] out = new double[]{sampleNumSimTotal, etaFrac, 1-serviceRateSAA, ObjSAA/100.0, objSAASim, SAAServiceSim, ObjScenario, ObjScenarioSim, 
//				ScenarioServiceSim, SDPSurvival, serviceSDP, SDPLbSurvival, serviceSDPLB, rollingObj, serviceRateRolling, rollingLength};
		
//		double[] out = new double[]{sampleNumSimTotal, etaFrac, 1-serviceRateSAA, 0, 0, 0, 0, 0, 
//				0, 0, 0, 0, 0, rollingObj, serviceRateRolling, rollingLength};
		
//		wr.writeToExcelAppend(out, fileName);	
//		System.out.println("");
//		
//		}
	}
}





