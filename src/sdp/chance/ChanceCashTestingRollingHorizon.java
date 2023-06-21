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


import milp.LostSaleChanceTesting;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.RiskRecursion;
import sdp.cash.RiskSimulation;
import sdp.cash.RiskState;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.CartesianProduct;
import sdp.sampling.Sampling;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.charts.SSJCategorySeriesCollection;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.rng.RandomStream;

/**
 * @author chen
 * @email: okchen321@163.com
 * @date: 2022 Mar 26, 15:14:50  
 * @desp: 
 * 
 * test the effects of sample number, service level, planning horizon length, distribution.
 * simulate 6 periods problem of SAA method may take more than 40 minutes.
 * 
 * increase the number of sample number will increase the performance of rolling horizon;
 * increase the rolling length may increase the performance;
 * increase the service rate will increase the feasibility.
 *
 */
public class ChanceCashTestingRollingHorizon {
	
	public static void main(String[] args) {
		WriteToExcelTxt wr = new WriteToExcelTxt();
		String fileName = "RollingHorizonTest6Periods-3.xls";
		String headString =  
				"demand mode" +
//						+ "\t" + "SDP survival" + "\t" + "SDP service" + "\t" + "SDP time" + 
//						"\t" + "INSDP survival" + "\t" + "InSDP service" + "\t" + "InSDP time" + 
						"\t"+ "rolling length" + "\t"  
						+ "serviceRate" + "\t" + "sample one period" +  "\t"+ "scenarioNumRolling" 
						+ "\t" + "iniCash" + "\t" + "price" + "\t" + 
						"variCost" + "\t" + "rolling time" + "\t" + "rolling obj" + "\t" + "rolling service rate";
		wr.writeToFile(fileName, headString);
		
		double[][] meanDemands = {{30,30,30,30,30,30,30, 30, 30, 30, 30, 30},
				{46,49,50,50,49,46,42,38,33,28,23,18},
				{11, 14, 18, 23, 28, 33, 38, 42, 46, 49, 50, 49},
				{47,30,13,6,13,30,47,54,47,30,13,6},
				{36,30,24,21,24,30,36,39,36,30,24,21},
				{63,27,10,24,1,23,33,35,67,7,14,41},
				{5,15,46,140,80,147,134,74,84,109,47,88},
				{14,24,71,118,49,86,152,117,226,208,78,59},
				{13,35,79,43,44,59,22,55,61,34,50,95},
				{15,56,19,84,136,67,67,155,87,164,19,67}
		};
		
//		double[][] meanDemands = {{30,30,30,30,30,30},
//				{50, 6,	38,	28,	23,	18},
//				{14,18,	23,	33,	42,	49},
//				{47,30,13,30,47,54},
//				{21,24,	39,	30,	24,	21},
//				{63,10,	4,	33,	67,	14},
//				{15,140,147,74,	109,88},
//				{14,71,	49,	152,78,	33},
//				{13,35,	79,	43,	44,	59},
//				{15,56,	19,	84,	136,67}
//		};
						
		double iniI = 0;
		double trunQuantile = 0.99;
		double[] serviceRate = new double[]{0.8}; // the higher value results in slower running speed. maximum negative possibility rate is 1 - serviceRate. 
		
		int T = meanDemands[0].length;
		int[] sampleNums = new int[T];
		double[] prices = new double[T];
		double[] variCostUnits = new double[T];
		double[] overheadCosts = new double[T];	
		double[] seasonalPrice = {14, 20, 18, 15};
		double[] seasonalVariCost = {10, 15, 16, 12};
		double sigmasCoe = 0.25; // for use in normal distribution
		
		
//		for (int t = 0; t < T; t++) {
//			int m = t % 4;
//			prices[t] = seasonalPrice[m];
//			variCostUnits[t] = seasonalVariCost[m];
//		}
		Arrays.fill(prices, 5);	
		Arrays.fill(variCostUnits, 1); 
		Arrays.fill(overheadCosts, 80); // overhead costs
		
		double holdCostUnit = 0;
		double salvageValueUnit = 0.5;
		double maxOrderQuantity = 1000; // maximum ordering quantity
		
		int[] sampleNumsRolling = new int[T];
		int sampleOnePeriod; 
		
		int rollingLength = 3; // rolling horizon length		
		int instanceNum = meanDemands.length;
		double iniCash = 130;
		int sampleNumRolling = 100; // number of scenarios for rolling SAA
		
		DecimalFormat df = new DecimalFormat("###, ###");
		NumberFormat nf = NumberFormat.getPercentInstance();
		nf.setMinimumFractionDigits(5);
		
		for (int kSampleNum = 0; kSampleNum < 3; kSampleNum ++)
			for (int kRollingLength = 1; kRollingLength < 2; kRollingLength++)
				for (int kService = 1; kService < 2; kService ++) {
					rollingLength = kRollingLength + 1;
					sampleOnePeriod = 10 + kSampleNum * 20;
					serviceRate[0] = 0.9 + kService * 0.05;	
					Arrays.fill(sampleNumsRolling,sampleOnePeriod);
					
					
		for (int m = 7; m < 8; m++) {
			for(int kk = 0; kk < 1; kk++) {	
					
			int instanceIndex = m;
			Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
						//.mapToObj(i -> new NormalDist(meanDemands[instanceIndex][i], sigmasCoe*meanDemands[instanceIndex][i])).toArray(Distribution[]::new); //can be changed to other distributions
						.mapToObj(i -> new PoissonDist(meanDemands[instanceIndex][i])).toArray(Distribution[]::new); // Poisson demand
			
			/**
		     * for simulation class
			 */
		    int stepSize = 1;
		    double fixOrderCost = 0;
		    double depositeRate = 0;
		    double minInventoryState = 0;
			double maxInventoryState = 500;
			double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
			double maxCashState = 3000;
			double discountFactor = 1;
			double[][][] pmf = new GetPmf(distributions, trunQuantile, stepSize).getpmf();	
			
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
			long currTime = System.currentTimeMillis();
			
			double SDPSurvival = 0;
			double serviceSDP = 0;
			double timeSDP = 0;
			double SDPLbSurvival = 0;
			double serviceSDPLb = 0;
			double timeSDPLb = 0;
			
//			if (kk == 0) {
//			/**
//			 * solve the problem by SDP when there is no joint chance constraint
//			 */	
//			// feasible actions
//			Function<RiskState, double[]> getFeasibleAction = s -> {
//				int t = s.getPeriod() - 1;
//				double maxQ = Math.min(s.iniCash/variCostUnits[t], maxOrderQuantity);
//				if (s.getBankruptBefore() == true)
//					maxQ = 0;
//				maxQ = Math.max(maxQ, 0);
//				return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
//			};	
//			
//			/*******************************************************************
//			 * Solve
//			 */
//			RiskRecursion recursion = new RiskRecursion(pmf, getFeasibleAction, stateTransition, immediateValue);
//			
//			
//			recursion.setTreeMapCacheAction();
//			SDPSurvival = recursion.getSurvProb(initialState);
//			System.out.println("**********************************************");
//			System.out.println("result of SDP is: ");
//			System.out.println("survival probability for this initial state is: " + nf.format(SDPSurvival));
//			double optQ = recursion.getAction(initialState);
//			System.out.println("optimal order quantity in the first priod is : " + optQ);
//			timeSDP = (System.currentTimeMillis() - currTime) / 1000;
//			System.out.println("running time is " + timeSDP + "s");
//			
//			/*******************************************************************
//			 * Simulating sdp results
//			 */	
//			int sampleNumSim = 1000;
//			RiskSimulation simulation = new RiskSimulation(distributions, sampleNumSim, recursion); // no need to add overheadCost in this class
//			double[] result2 = simulation.simulateLostSale(initialState, immediateValue);
//			System.out.println("final simulated survival probability in " + df.format(sampleNumSim) + " samples is: " + nf.format(result2[0]));
//			System.out.println("final simulated lost sale rate " + " is: " + nf.format(result2[1])); 
//			serviceSDP = 1 - result2[1];
//			System.out.println("final simulated service sale rate " + " is: " + nf.format(serviceSDP)); 
//			
//			/*******************************************************************
//			 * solve the problem by SDP when there is individual chance constraint approximation
//			 */			
//			// feasible actions 2
//			Function<RiskState, double[]> getFeasibleAction2 = s -> {
//				int t = s.getPeriod() - 1;
//				double thisPeriodServRate = 1- (1 - serviceRate[0]) / T;
//				double minQ = Math.ceil(distributions[t].inverseF(thisPeriodServRate)); // minimum ordering quantity in each period 
//				double maxQ = Math.min(s.iniCash/variCostUnits[t], maxOrderQuantity);
//				if (s.getBankruptBefore() == true)
//					maxQ = 0;
//				if (maxQ < minQ) {
//					maxQ = 0;
//					minQ = 0;
//				}
//				maxQ = Math.max(maxQ, 0);
//				return DoubleStream.iterate(minQ, i -> i + stepSize).limit((int) maxQ + 1).toArray();
//			};
//			
//			/*******************************************************************
//			 * Solve
//			 */
//			recursion = new RiskRecursion(pmf, getFeasibleAction2, stateTransition, immediateValue);
//			period = 1;		
//			initialState = new RiskState(period, iniI, iniCash, false);
//			currTime = System.currentTimeMillis();
//			recursion.setTreeMapCacheAction();
//			
//			SDPLbSurvival= recursion.getSurvProb(initialState);
//			System.out.println("**********************************************");
//			System.out.println("result of SDP with service rate constraint is: ");
//			System.out.println("survival probability for this initial state is: " + nf.format(SDPLbSurvival));
//			System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
//			timeSDPLb = (System.currentTimeMillis() - currTime) / 1000;
//			System.out.println("running time is " + timeSDPLb + "s");
//			
//			/*******************************************************************
//			 * Simulating sdp results
//			 */		
//			sampleNumSim = 10000;
//			simulation = new RiskSimulation(distributions, sampleNumSim, recursion); // no need to add overheadCost in this class
//			result2 = simulation.simulateLostSale(initialState, immediateValue);
//			System.out.println("final simulated survival probability in " + df.format(sampleNumSim) + " samples is: " + nf.format(result2[0]));
//			System.out.println("final simulated lost sale rate " + " is: " + nf.format(result2[1]));
//			serviceSDPLb = 1 - result2[1];
//			System.out.println("final simulated service rate " + " is: " + nf.format(serviceSDPLb));
//			}
			
			/**
			 * solve the problem by rolling horizon of SAA
			 * 
			 */
			RiskSimulation simulation1 = new RiskSimulation(distributions, sampleNumRolling, immediateValue, stateTransition); // no need to add overheadCost in this class    		
		    currTime = System.currentTimeMillis();
		    System.out.println("**********************************************");
		    System.out.println("after rolling horizon for length " + rollingLength +", result of SAA: ");
		    double[] resultRolling; 
		    Sampling sampling = new Sampling();
		    
		    double meanDemandSum = Arrays.stream(meanDemands[m]).sum();
			double rollingDemandSum = Arrays.stream(meanDemands[m]).limit(rollingLength).sum();
			double portion = rollingDemandSum / meanDemandSum;
			double rollingServiceRate = Math.pow(serviceRate[0], portion);
				
			
		    double[][] scenariosRolling = sampling.generateLHSamples(distributions, sampleNumsRolling); //include scenarios in the whole planning horizon
		    resultRolling = simulation1.rollingHoirzon(rollingLength, initialState, rollingServiceRate, sampleNumsRolling, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenariosRolling, sampleNumRolling);
		    double time1 = (System.currentTimeMillis() - currTime) / 1000.00;	     
		    System.out.println("running time is " + time1 + "s");
		    System.out.println("each period sample number is " + sampleOnePeriod);
		    System.out.println("final simulated survival probability of rolling  SAA in " + df.format(sampleNumRolling) + "samples is: " + nf.format(resultRolling[0]));
		    double sigma2 = Math.sqrt(resultRolling[1]*(1 - resultRolling[1])/sampleNumRolling);
			double error2  = 1.96*sigma2;
			double serviceRate2 = 1 - resultRolling[1];
			System.out.println("the service rate for SAA rolling horizon is " + nf.format(serviceRate2)  + " with error " + nf.format(error2));
			
			/*******************************************************************
			 * 
			 * output results to excel
			 */
			System.out.println("");
			// SDPSurvival, serviceSDP, timeSDP, SDPLbSurvival, serviceSDPLb, timeSDPLb,
			double[] out = new double[]{m, rollingLength, serviceRate[0], sampleOnePeriod, sampleNumRolling, iniCash, prices[0], variCostUnits[0], time1, resultRolling[0], serviceRate2};
			wr.writeToExcelAppend(out, fileName);
			}
		}			
		}
	}
}