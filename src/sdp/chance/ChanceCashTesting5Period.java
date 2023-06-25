package sdp.chance;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.util.function.IntToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;


import milp.LostSaleChance;
import milp.LostSaleChanceTesting;
import milp.testPiecewise;
import sdp.cash.CashRecursion;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.RiskRecursion;
import sdp.cash.RiskSimulation;
import sdp.cash.RiskState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.CartesianProduct;
import sdp.sampling.Sampling;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.charts.SSJCategorySeriesCollection;
import umontreal.ssj.probdist.BernoulliDist;
import umontreal.ssj.probdist.BinomialDist;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.probdistmulti.BiNormalDist;
import umontreal.ssj.rng.RandomStream;

/**
 * @author chen
 * @email: okchen321@163.com
 * @date: 2022 Mar 26, 15:14:50  
 * @desp: 
 * 
 * test the effects of sample number, service level, planning horizon length, fluctuation, holding cost, initial cash/ margin.
 * simulate 6 periods problem of SAA method may take more than 40 minutes.
 * 
 *
 */
public class ChanceCashTesting5Period {
	
	public static void main(String[] args) {
		WriteToExcelTxt wr = new WriteToExcelTxt();
		String fileName = "RollingTest5Periods.xls";
		String headString =  
				"demand mode" + "\t" + "serviceRate" + "\t" + "sample number" + "\t" + "iniCash" + "\t" + "price" + "\t" +  
		         "overheadCost" + "\t" + "SDPObj" + "\t" + "SDPService" + "\t" + "timeSDP" + "\t" + "SDPLbObj" + "\t" 
					 + "SDPLbServie" + "\t" + "timeSDPLb" + "\t" + "RollingObj" + "\t" + "RollingService" +  "\t" + "RollingTime"
					 + "\t" + "rollingLength" + "\t" + "Q1SDP" + "\t" + "Q1SDPLb" + "\t" + "Q1Rolling" + "\t";
		wr.writeToFile(fileName, headString);
		
		
		double[][] meanDemands = {{30,30,30,30,30},
				{50, 46, 38, 28, 14},
				{14,23,33,46,50},
				{47,30,6,30,54},
				{9,30,44,30,8 },
				{63,27,10,24,1},
				{25, 46, 140, 80, 147},
				{14,24,71,118,49},
				{13,35,79,43,44},
				{15,56,19,84,136} };
		
		double[] iniCashT = {80, 100, 120};
		double[] pricesT = {4, 5, 6};
		double[] overheadCostsT = {60, 80, 100};
				
		double iniI = 0;
		double trunQuantile = 0.999;
		double serviceRate = 0.95; // the higher value results in slower running speed. maximum negative possibility rate is 1 - serviceRate. 
		
		int T = meanDemands[0].length;
		int[] sampleNums = new int[T];
		double[] prices = new double[T];
		double[] variCostUnits = new double[T];
		double[] overheadCosts = new double[T];	
		double[] seasonalPrice = {14, 20, 18, 15};
		double[] seasonalVariCost = {10, 15, 16, 12};
		double holdCostUnit = 0;
		
		double sigmasCoe = 0.25; // for use in normal distribution
		int sampleNum = 100;  // sample number in rolling horizon simulation
		
		for (int t = 0; t < T; t++) {
			int m = t % 4;
			prices[t] = seasonalPrice[m];
			variCostUnits[t] = seasonalVariCost[m];
		}
		Arrays.fill(prices, 16);	//
		Arrays.fill(variCostUnits, 10); //
		Arrays.fill(overheadCosts, 100); // overhead costs
		
		double salvageValueUnit = 0.5;
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash
				
		
		int instanceNum = meanDemands.length;
		int runNum = 1; // run number is 10 for computing upper bounds
		
		for (int iCash = 0; iCash < 3; iCash++)
			for (int iPrice = 0; iPrice < 3; iPrice++)
				for (int iOverh = 0; iOverh < 3; iOverh++) {
					double iniCash = iniCashT[iCash];
					Arrays.fill(prices, pricesT[iPrice]);
					Arrays.fill(overheadCosts, overheadCostsT[iOverh]);
					Arrays.fill(variCostUnits, 1);
					
					int sampleNumPeriod = 300; // number of samples in each period
					int rollingLength = 1; // rolling horizon length
		
		for(int kk = 0; kk < 1; kk++) {
		for (int m = 0; m < instanceNum; m++) {
			
			double time1;
			LostSaleChance model;
			double[] resultSAA = new double[3];			
			double positiveScenario;
			double survivalProb;
			double lostRate;
			NumberFormat nf = NumberFormat.getPercentInstance();
			nf.setMinimumFractionDigits(5);		
			DecimalFormat df = new DecimalFormat("###, ###");
			
			double[] sAA10Run = new double[runNum];
			double[] sAA10Time = new double[runNum];
			double[] Q10Run = new double[runNum];
			
			Arrays.fill(sampleNums, sampleNumPeriod);	
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			
			int instanceIndex = m;
			Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
						//.mapToObj(i -> new NormalDist(meanDemands[instanceIndex][i], sigmasCoe*meanDemands[instanceIndex][i])).toArray(Distribution[]::new); //can be changed to other distributions
						.mapToObj(i -> new PoissonDist(meanDemands[instanceIndex][i])).toArray(Distribution[]::new); // Poisson demand
		
			double meanDemandSum = Arrays.stream(meanDemands[m]).sum();
			double rollingDemandSum = Arrays.stream(meanDemands[m]).limit(rollingLength).sum();
			double portion = rollingDemandSum / meanDemandSum;
			double rollingServiceRate = Math.pow(serviceRate, portion);
			Sampling sampling = new Sampling();

			for (int k = 0; k < runNum; k++) {		
				
		    int stepSize = 1;
		    double fixOrderCost = 0;
		    double depositeRate = 0;
		    double minInventoryState = 0;
			double maxInventoryState = 800;
			double minCashState = -1000; // can affect results, should be smaller than minus fixedOrderCost
			double maxCashState = 2000;
			double discountFactor = 1;
			double[][][] pmf = new GetPmf(distributions, trunQuantile, stepSize).getpmf();	
			
			// feasible actions
			Function<RiskState, double[]> getFeasibleAction = s -> {
				int t = s.getPeriod() - 1;
				double maxQ = Math.min(s.iniCash/variCostUnits[t], maxOrderQuantity);
				if (s.getBankruptBefore() == true)
					maxQ = 0;
				maxQ = Math.max(maxQ, 0);
				return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
			};
			
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
						  
			int period = 1;	
			RiskState initialState = new RiskState(period, iniI, iniCash, false);			    			    			    	    
			RiskSimulation simulation1 = new RiskSimulation(distributions, sampleNum, immediateValue, stateTransition); // no need to add overheadCost in this class
		    	
		    
		    long currTime = System.currentTimeMillis();
			
			/*******************************************************************
			 * Solve SDP
			 */
			RiskRecursion recursion = new RiskRecursion(pmf, getFeasibleAction, stateTransition, immediateValue);
			recursion.setTreeMapCacheAction();
			double SDPObj = recursion.getSurvProb(initialState);
			System.out.println("survival probability for this initial state is: " + SDPObj);
			double optQSDP = recursion.getAction(initialState);
			System.out.println("optimal order quantity in the first priod is : " + optQSDP);
			double timeSDP = (System.currentTimeMillis() - currTime) / 1000;
			System.out.println("running time is " + timeSDP + "s");
			
			System.out.println();
			
			/*******************************************************************
			 * Simulate the result
			 */
			
			int sampleNumSim = 1000;
			RiskSimulation simulation = new RiskSimulation(distributions, sampleNumSim, recursion); // no need to add overheadCost in this class
			double[] result = simulation.simulateLostSale(initialState, immediateValue);
			DecimalFormat df2 = new DecimalFormat("###, ###");
			System.out.println("\nfinal simulated survival probability in " + df2.format(sampleNumSim) + " samples is: " + nf.format(result[0]));
			double SDPService = 1 - result[1];
			System.out.println("final simulated service rate " + " is: " + nf.format(SDPService));
			System.out.println("************************************************************");
					
			/*******************************************************************
			 * solve the problem by SDP when there is individual chance constraint approximation,
			 * to get a lower bound
			 */			
			Function<RiskState, double[]> getFeasibleAction2 = s -> {
				int t = s.getPeriod() - 1;
				double thisPeriodServRate = 1 - (1 - serviceRate) / T;
				double minQ = Math.ceil(distributions[t].inverseF(thisPeriodServRate)); // minimum ordering quantity in each period 
				double maxQ = Math.min(s.iniCash/variCostUnits[t], maxOrderQuantity);
				if (s.getBankruptBefore() == true)
					maxQ = 0;
				if (maxQ < minQ) {
					maxQ = 0;
					minQ = 0;
				}
				maxQ = Math.max(maxQ, 0);
				return DoubleStream.iterate(minQ, i -> i + stepSize).limit((int) maxQ + 1).toArray();
			};
			
			/*******************************************************************
			 * Solve
			 */		
			recursion = new RiskRecursion(pmf, getFeasibleAction2, stateTransition, immediateValue);
			period = 1;		
			initialState = new RiskState(period, iniI, iniCash, false);
			currTime = System.currentTimeMillis();
			recursion.setTreeMapCacheAction();
			double SDPLbObj;
			SDPLbObj = recursion.getSurvProb(initialState);
			System.out.println("result of SDP with service rate constraint is: ");
			System.out.println("survival probability for this initial state is: " + nf.format(SDPLbObj));
			double optQSDPLb = recursion.getAction(initialState);
			System.out.println("optimal order quantity in the first priod is : " + optQSDPLb);
			double timeSDPLb = (System.currentTimeMillis() - currTime) / 1000;
			System.out.println("running time is " + timeSDPLb + "s");
			
			/*******************************************************************
			 * Simulating sdp results
			 */		
			
			sampleNum = 1000;
			simulation = new RiskSimulation(distributions, sampleNum, recursion); // no need to add overheadCost in this class
			result = simulation.simulateLostSale(initialState, immediateValue);
			System.out.println("final simulated survival probability in " + df.format(sampleNum) + " samples is: " + nf.format(result[0]));
//			System.out.println("final simulated lost sale rate " + " is: " + nf.format(result[1]));
			double serviceSDPLB = 1 - result[1];
			System.out.println("final simulated service rate " + " is: " + nf.format(serviceSDPLB));
			
			/*******************************************************************
			 * rolling horizon approach
			 * 
			 */
			int[] sampleNumSims = new int[T]; // sample number in each period, the number of samples in the first period can have big influence		
			int[] sampleNumSimsRolling = new int[T];
			Arrays.fill(sampleNumSimsRolling, sampleNumPeriod);
		    currTime = System.currentTimeMillis();
		    System.out.println("**********************************************");
		    
		    // generate samples
		 	double[][] scenarios = sampling.generateLHSamples(distributions, sampleNumSims);
		    scenarios = sampling.generateLHSamples(distributions, sampleNumSimsRolling);
		    
		   
		    
		    double[] resultSim; 	
		    int sampleNumRolling = 1000;
		    resultSim = simulation1.rollingHoirzon(rollingLength, initialState, rollingServiceRate, sampleNumSimsRolling, prices, variCostUnits, 
		    					overheadCosts, salvageValueUnit, holdCostUnit, scenarios, sampleNumRolling);
		    double timeRolling = (System.currentTimeMillis() - currTime) / 1000.00;	    
		    System.out.println("after rolling horizon for length " + rollingLength +", " + "total horizon length is " + T + ", result is: ");
		    System.out.println("running time is " + timeRolling + "s");
		    double rollingObj = resultSim[0];
		    System.out.println("final simulated survival probability of rolling SAA in " + sampleNumRolling + " samples is: " + nf.format(rollingObj));
		    double sigma2 = Math.sqrt(resultSim[1]*(1 - resultSim[1])/sampleNum);
			double error2  = 1.96*sigma2;
			double serviceRateRolling = 1 - resultSim[1];
			double optQ1Rolling = resultSim[2];
			System.out.println("the service rate for simulated SAA rolling horizon is " + nf.format(serviceRateRolling) + ", with error " + nf.format(error2));
			System.out.println("the optimal ordering Q in the 1st period of the SAA rolling horizon is " + optQ1Rolling);
			System.out.println("**********************************************");
			System.out.println("**********************************************");
			
			/*******************************************************************
			 * 
			 * output results to excel
			 */
			double[] out = new double[]{m, serviceRate, sampleNumPeriod, iniCash, prices[0], overheadCosts[0], SDPObj, SDPService, timeSDP, SDPLbObj, serviceSDPLB, timeSDPLb,
					rollingObj, serviceRateRolling, timeRolling, rollingLength, optQSDP, optQSDPLb, optQ1Rolling};
			
			wr.writeToExcelAppend(out, fileName);
			}
		}
		}
			    
		}
	}
}