package sdp.chance;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.IntToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


import milp.LostSaleChanceTesting;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.CartesianProduct;
import sdp.sampling.Sampling;
import sdp.write.WriteToExcel;
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
 * need seasonal price and vari cost
 *
 */
public class ChanceCashTesting {
	
	public static void main(String[] args) {
		WriteToExcel wr = new WriteToExcel();
		String fileName = "JointChanceSAA.xls";
		String headString =  
				"demand mode" + "\t" + "SAA obj" + "\t" + "time" + "t" + "sim SAA obj" + "\t" + "sim obj error" + "\t" +
		         "sim SAA service rate" + "\t" + "sim service error" + "sim time" + "\t" + "extend SAA obj" + "\t" + "time" + "t" + "sim extend SAA obj" + "\t" 
					+ "t" + "obj error" + "t" + "sim extend SAA service" + "\t" + "service error" + "t" + "sim time";
		
		double[][] meanDemands = {{30,30,30,30,30,30,30, 30, 30, 30, 30, 30},
				{46,49,50,50,49,46,42,38,33,28,23,18},
				{8, 11, 14, 18, 23, 28, 33, 38, 42, 46, 49, 50},
				{47,30,13,6,13,30,47,54,47,30,13,6},
				{36,30,24,21,24,30,36,39,36,30,24,21},
				{63,27,10,24,1,23,33,35,67,7,14,41},
				{1,15,46,140,80,147,134,74,84,109,47,88},
				{14,24,71,118,49,86,152,117,226,208,78,59},
				{13,35,79,43,44,59,22,55,61,34,50,95},
				{15,56,19,84,136,67,67,155,87,164,19,67}
		};
		
		double iniCash = 40;
		double iniI = 0;
		double trunQuantile = 0.999;
		double serviceRate = 0.9; // the higher value results in slower running speed. maximum negative possibility rate is 1 - serviceRate. 
		
		int T = 6; //meanDemands[0].length;
		int[] sampleNums = new int[T];
		double[] prices = new double[T];
		double[] variCostUnits = new double[T];
		double[] overheadCosts = new double[T];	
		double[] seasonalPrice = {14, 20, 18, 15};
		double[] seasonalVariCost = {10, 15, 16, 12};
		double sigmasCoe = 0.25; // for use in normal distribution
		int sampleTotalNumber = 3000;
		int sampleNum = 100;  // sample number in simulation
		
		int sampleNumPeriod = 5; // number of samples in each period
		Arrays.fill(sampleNums, sampleNumPeriod);
		for (int t = 0; t < T; t++) {
			int m = t % 4;
			prices[t] = seasonalPrice[m];
			variCostUnits[t] = seasonalVariCost[m];
		}
		Arrays.fill(prices, 20);	
		Arrays.fill(variCostUnits, 10);
		Arrays.fill(overheadCosts, 100); // overhead costs
		double holdCostUnit = 1;
		double salvageValueUnit = 5;
		double maxOrderQuantity = 300; // maximum ordering quantity when having enough cash
		
		int[] sampleNumsRolling = new int[T];
		Arrays.fill(sampleNumsRolling,5);
		int rollingLength = 4; // rolling horizon length		
		
		int instanceNum = meanDemands.length;
		int runNum = 10;
		
		for (int m = 0; m < instanceNum; m++) {
			for (int k = 0; k < runNum; k++) {					
				int instanceIndex = m;
				Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
						//.mapToObj(i -> new NormalDist(meanDemand[i], sigmasCoe*meanDemand[i])).toArray(Distribution[]::new); //can be changed to other distributions
						.mapToObj(i -> new PoissonDist(meanDemands[instanceIndex][i])).toArray(Distribution[]::new); // Poisson demand
		
				double meanDemandSum = Arrays.stream(meanDemands[m]).sum();
				double rollingDemandSum = Arrays.stream(meanDemands[m]).limit(rollingLength).sum();
				double portion = rollingDemandSum / meanDemandSum;
				double rollingServiceRate = Math.pow(serviceRate, portion);
				
				// generate scenarios, samples in each period form a scenario tree
				Sampling sampling = new Sampling();
				double[][] scenarios = sampling.generateLHSamples(distributions, sampleNums);
		
				// sample without replacement
				Random rd = new Random();
				double[][] demandSamples = new double[sampleTotalNumber][T];
				for (int i = 0; i < sampleTotalNumber; i++) {
					for (int t = 0; t < T; t++) {
						int index = rd.nextInt(sampleNumPeriod);
						demandSamples[i][t] = scenarios[t][index];
					}
				}
				
				/**
				 * solve the problem by SAA
				 */				
				int negativeScenarioNumRequire = (int) (sampleTotalNumber * (1 - serviceRate));
				LostSaleChanceTesting model = new LostSaleChanceTesting(distributions, demandSamples, iniCash, iniI, prices, variCostUnits, 
									salvageValueUnit, holdCostUnit, overheadCosts, serviceRate);
				
				long currTime = System.currentTimeMillis();
				double[] resultSAA;
				double time1;
				double positiveScenario;
				double survivalProb;
				double lostRate;
				NumberFormat nf = NumberFormat.getPercentInstance();
				nf.setMinimumFractionDigits(5);		
				DecimalFormat df = new DecimalFormat("###, ###");
				
				resultSAA = model.solveMaxSurvivalTesting();				
				time1 = (System.currentTimeMillis() - currTime) / 1000.00;
			    currTime = System.currentTimeMillis();	    
			    System.out.println("**********************************************");
			    System.out.println("result of SAA-scenario tree before sorting scenarios: ");
			    System.out.println("running time is " + time1 + "s");	
			    System.out.printf("first stage decison Q is: %.2f\n", resultSAA[0]);
			    positiveScenario = resultSAA[1];
			    System.out.printf("Objective value is: %.0f in %d scenarios\n", resultSAA[1], sampleTotalNumber);
			    survivalProb = 100 * resultSAA[1] / sampleTotalNumber;
			    System.out.printf("Survival probability is: %.5f%%\n", survivalProb);
			    System.out.println("lost sale scenario number in the solution is : " + resultSAA[2]);
			    System.out.println("maximum lost sale scenario number allowed is: " + negativeScenarioNumRequire);
			    lostRate = resultSAA[2] / (double) sampleTotalNumber;
			    System.out.println("lost sale rate of SAA is: " + nf.format(lostRate));
			    System.out.println("lost sale max required rate is: " + nf.format(1 - serviceRate));
			    System.out.println();
			    
				/**
				 * solve the problem by extended formulation of SAA,
				 * sort scenarios in the whole planning horizon
				 * 
				 */
				currTime = System.currentTimeMillis();
				double[] resultExtendSAA;
				resultExtendSAA = model.solveSortWholeTesting();	// same result with soveSort or solveSort2, but less computational time
				                                 				
				time1 = (System.currentTimeMillis() - currTime) / 1000.00;
			    currTime = System.currentTimeMillis();	    
			    System.out.println("**********************************************");
			    System.out.println("after sorting scenarios in the whole planning horizon, result of SAA-scenario tree: ");
			    System.out.println("running time is " + time1 + "s");	
			    System.out.printf("first stage decison Q is: %.2f\n", resultExtendSAA[0]);
			    positiveScenario = resultExtendSAA[1];
			    System.out.printf("Objective value is: %.0f in %d scenarios\n", positiveScenario, sampleTotalNumber);
			    survivalProb = 100 * resultExtendSAA[1] / sampleTotalNumber;
			    System.out.printf("Survival probability is: %.5f%%\n", survivalProb);
			    System.out.println("lost sale scenario number in the solution is : " + resultExtendSAA[2]);
			    System.out.println("maximum lost sale scenario number allowed is: " + negativeScenarioNumRequire);
			    lostRate = resultExtendSAA[2] / (double) sampleTotalNumber;
			    System.out.println("lost sale rate of SAA is: " + nf.format(lostRate));
			    System.out.println("lost sale max required rate is: " + nf.format(1 - serviceRate));
			    System.out.println();
			    
			    /**
			     * for simulation class
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
			    
				/**
				 * solve the problem by rolling horizon of extended SAA
				 * 
				 */
				sampleNum = 100; // number of scenarios for rolling SAA
			    currTime = System.currentTimeMillis();
			    System.out.println("**********************************************");
			    System.out.println("after rolling horizon for length " + rollingLength +", result of SAA-scenario tree: ");
			    double[] resultRolling; 
			    double[][] scenariosRolling = sampling.generateLHSamples(distributions, sampleNums); //include scenarios in the whole planning horizon
			    resultRolling = simulation1.rollingHoirzonFurtherExtendSAA(rollingLength, initialState, rollingServiceRate, sampleNumsRolling, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenariosRolling, sampleNum);
			    time1 = (System.currentTimeMillis() - currTime) / 1000.00;	     
			    System.out.println("running time is " + time1 + "s");
			    System.out.println("final simulated survival probability of rolling further SAA in " + df.format(sampleNum) + " samples is: " + nf.format(resultRolling[0]));
		    		
			    /**
			     * Simulate the restult of SAA
				 */			    
			    double error;
			    double thisServiceRate;		
			    double[] result1; 
			    result1 = simulation1.simulateSAATesting(initialState, resultSAA[0], serviceRate, scenarios, sampleNums, demandSamples, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, sampleNum);
				System.out.println("final simulated survival probability of SAA in " + df.format(sampleNum) + " samples is: " + nf.format(result1[0]));
				error  = 1.96 * Math.sqrt(result1[1]*(1 - result1[1]) / sampleNum);
				thisServiceRate = 1-result1[1];
				System.out.println("final simulated service sale rate of SAA " + " is: " + nf.format(thisServiceRate) + " with error " + nf.format(error)); 
				System.out.println();
				
				/**
			     * Simulate the result of extended SAA
				 */	
			    result1 = simulation1.simulateExtendSAAWholeTesting(initialState, resultExtendSAA[0], serviceRate, demandSamples, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, sampleNum);
				System.out.println("final simulated survival probability of extended SAA(sort whole planning horizon) in " + df.format(sampleNum) + " samples is: " + nf.format(result1[0]));
				error  = 1.96 * Math.sqrt(result1[1]*(1 - result1[1]) / sampleNum);
				thisServiceRate = 1-result1[1];
				System.out.println("final simulated service rate of extended SAA(sort whole planning horizon) " + " is: " + nf.format(thisServiceRate) + " with error " + nf.format(error));

			}
		}
		System.out.println();
		
		
	}

}
