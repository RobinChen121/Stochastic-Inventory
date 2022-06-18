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
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.CartesianProduct;
import sdp.sampling.Sampling;
import sdp.write.WriteToExcel;
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
		WriteToExcel wr = new WriteToExcel();
		String fileName = "JointChanceSAA5Periods.xls";
		String headString =  
				"demand mode" + "\t" + "serviceRate" + "\t" + "scenario number" + "\t" + "iniCash" + "\t" + "price" + "\t" + "variCost" + "\t" + "SAA obj" + "\t" + "time" + "\t" + "sim SAA obj" + "\t" +
		         "sim SAA service rate" + "\t" + "sim saa time" + "\t" + "extend SAA obj" + "\t" + "time" + "\t" + "sim extend SAA obj" + "\t" 
					 + "sim extend SAA service" + "\t" + "sim extend time" + "\t" + "upper bound" + "\t" + "lower bound" +  "\t" + "simSampleNum"
					 + "\t" + "sigmasCoe" + "\t" +"holdingCost";
		wr.writeToFile(fileName, headString);
		
		
		double[][] meanDemands = {{30,30,30,30,30},
				{50, 46, 38, 28, 14},
				{14,23,33,46,50},
				{47,30,6,30,54},
				{9,30,44,30,8 },
				{63,27,10,24,1},
				{15, 46, 140, 80, 147},
				{14,24,71,118,49},
				{13,35,79,43,44},
				{15,56,19,84,136} };
		
		double iniCash = 40;
		double iniI = 0;
		double trunQuantile = 0.999;
		double serviceRate = 0.6; // the higher value results in slower running speed. maximum negative possibility rate is 1 - serviceRate. 
		
		int T = meanDemands[0].length;
		int[] sampleNums = new int[T];
		double[] prices = new double[T];
		double[] variCostUnits = new double[T];
		double[] overheadCosts = new double[T];	
		double[] seasonalPrice = {14, 20, 18, 15};
		double[] seasonalVariCost = {10, 15, 16, 12};
		double sigmasCoe = 0.25; // for use in normal distribution
		int sampleNum = 100;  // sample number in simulation
		
		for (int t = 0; t < T; t++) {
			int m = t % 4;
			prices[t] = seasonalPrice[m];
			variCostUnits[t] = seasonalVariCost[m];
		}
		Arrays.fill(prices, 16);	//
		Arrays.fill(variCostUnits, 10); //
		Arrays.fill(overheadCosts, 100); // overhead costs
		double holdCostUnit = 0;
		double salvageValueUnit = 2;
		double maxOrderQuantity = 300; // maximum ordering quantity when having enough cash
		
		int[] sampleNumsRolling = new int[T];
		Arrays.fill(sampleNumsRolling,5);
		int rollingLength = 3; // rolling horizon length		
		
		int instanceNum = meanDemands.length;
		int runNum = 1; // run number is 10 for computing upper bounds
		
		for(int kk = 0; kk < 10; kk++) {
		for (int m = 0; m < instanceNum; m++) {
			long currTime = System.currentTimeMillis();
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
			int sampleNumPeriod = 5; // number of samples in each period
			Arrays.fill(sampleNums, sampleNumPeriod);	
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			
			int instanceIndex = m;
			Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
						.mapToObj(i -> new NormalDist(meanDemands[instanceIndex][i], sigmasCoe*meanDemands[instanceIndex][i])).toArray(Distribution[]::new); //can be changed to other distributions
						//.mapToObj(i -> new PoissonDist(meanDemands[instanceIndex][i])).toArray(Distribution[]::new); // Poisson demand
		
			double meanDemandSum = Arrays.stream(meanDemands[m]).sum();
			double rollingDemandSum = Arrays.stream(meanDemands[m]).limit(rollingLength).sum();
			double portion = rollingDemandSum / meanDemandSum;
			double rollingServiceRate = Math.pow(serviceRate, portion);
			Sampling sampling = new Sampling();
			double[][] scenarios = new double[T][];
			for (int k = 0; k < runNum; k++) {					
				// generate scenarios, samples in each period form a scenario tree				
				scenarios = sampling.generateLHSamples(distributions, sampleNums);
				
				/**
				 * solve the problem by SAA
				 */				
				model = new LostSaleChance(distributions, sampleNums, iniCash, iniI, prices, variCostUnits, 
						salvageValueUnit, holdCostUnit, overheadCosts, serviceRate, scenarios);	
		        resultSAA = model.solveMaxSurvival();				
				time1 = (System.currentTimeMillis() - currTime) / 1000.00;
			    currTime = System.currentTimeMillis();	    
			    System.out.println("**********************************************");
			    System.out.println("result of SAA-scenario tree before sorting scenarios: ");
			    System.out.println("running time is " + time1 + "s");	
			    System.out.printf("first stage decison Q is: %.2f\n", resultSAA[0]);
			    positiveScenario = resultSAA[1];
			    System.out.printf("Objective value is: %.0f in %d scenarios\n", resultSAA[1], sampleNumTotal);
			    survivalProb = 100 * resultSAA[1] / sampleNumTotal;
			    System.out.printf("Survival probability is: %.5f%%\n", survivalProb);
			    System.out.println("lost sale scenario number in the solution is : " + resultSAA[2]);
			    System.out.println("maximum lost sale scenario number allowed is: " + negativeScenarioNumRequire);
			    lostRate = resultSAA[2] / (double) sampleNumTotal;
			    System.out.println("lost sale rate of SAA is: " + nf.format(lostRate));
			    System.out.println("lost sale max required rate is: " + nf.format(1 - serviceRate));
			    System.out.println();
			    sAA10Run[k] = survivalProb / 100.0;
			    sAA10Time[k] = time1;
			    Q10Run[k] = resultSAA[0];
			}			    
			Arrays.sort(sAA10Run);
			double thetaN = BinomialDist.cdf(sampleNumTotal, 1-serviceRate, negativeScenarioNumRequire);			
			double temp = BinomialDist.cdf(1024, 0.2, 204);
			double tau = 0.05; 
			int Lindex = 0;
			for (int i = 0; i < runNum; i++) {
				double prob = BinomialDist.cdf(runNum, thetaN, i);
				if ( prob > tau) {
					Lindex = i;
					break;
				}
			}
			double staUpperBound = sAA10Run[runNum - Lindex - 1];
			double saaObj = Arrays.stream(sAA10Run).average().getAsDouble();
			double saaTime = Arrays.stream(sAA10Time).average().getAsDouble();
					
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
		    
		    double error;
		    double thisServiceRate;		
		    double[] result1; 
		    
		    /**
		     * Simulate the restult of SAA
			 */			    			    
		    currTime = System.currentTimeMillis();
		    double Q = Arrays.stream(Q10Run).max().getAsDouble();
		    result1 = simulation1.simulateSAA(initialState, Q, serviceRate, sampleNums, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenarios, sampleNum);
		    time1 = (System.currentTimeMillis() - currTime) / 1000.00;  
		    System.out.println("running time is " + time1 + "s");
		    System.out.println("final simulated survival probability of SAA in " + df.format(sampleNum) + " samples is: " + nf.format(result1[0]));
			error  = 1.96 * Math.sqrt(result1[1]*(1 - result1[1]) / sampleNum);
			thisServiceRate = 1-result1[1];
			System.out.println("final simulated service sale rate of SAA " + " is: " + nf.format(thisServiceRate) + " with error " + nf.format(error)); 
			
			double saaSimObj = result1[0];
			double saaSimError = error;
			double saaSimTime = time1;
			double saaSimServiceRate = thisServiceRate;
			
			
			/**
			* solve the problem by extended formulation of SAA,
			* sort scenarios in the whole planning horizon
			* 
			*/
			currTime = System.currentTimeMillis();
			double[] resultExtendSAA;
			scenarios = sampling.generateLHSamples(distributions, sampleNums);
			model = new LostSaleChance(distributions, sampleNums, iniCash, iniI, prices, variCostUnits, 
					salvageValueUnit, holdCostUnit, overheadCosts, serviceRate, scenarios);	
			resultExtendSAA = model.solveSortWhole();	// same result with soveSort or solveSort2, but less computational time
				                                 				
			time1 = (System.currentTimeMillis() - currTime) / 1000.00;
		    currTime = System.currentTimeMillis();	    
		    System.out.println("**********************************************");
		    System.out.println("after sorting scenarios in the whole planning horizon, result of SAA-scenario tree: ");
		    System.out.println("running time is " + time1 + "s");	
		    System.out.printf("first stage decison Q is: %.2f\n", resultExtendSAA[0]);
		    positiveScenario = resultExtendSAA[1];
		    System.out.printf("Objective value is: %.0f in %d scenarios\n", positiveScenario, sampleNumTotal);
		    survivalProb = 100 * resultExtendSAA[1] / sampleNumTotal;
		    System.out.printf("Survival probability is: %.5f%%\n", survivalProb);
		    System.out.println("lost sale scenario number in the solution is : " + resultExtendSAA[2]);
		    System.out.println("maximum lost sale scenario number allowed is: " + negativeScenarioNumRequire);
		    lostRate = resultExtendSAA[2] / (double) sampleNumTotal;
		    System.out.println("lost sale rate of SAA is: " + nf.format(lostRate));
		    System.out.println("lost sale max required rate is: " + nf.format(1 - serviceRate));
		    System.out.println();
		    
		    double extendObj = survivalProb / 100.0;
		    double extendTime = time1;
			
		    /**
		     * Simulate the result of extended SAA
			 */	
			currTime = System.currentTimeMillis();
			System.out.println("**********************************************");
			result1 = simulation1.simulateExtendSAAWhole(initialState, resultExtendSAA[0], serviceRate, sampleNums, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenarios, sampleNum);
		    time1 = (System.currentTimeMillis() - currTime) / 1000.00;
		    System.out.println("running time is " + time1 + "s");
		    System.out.println("final simulated survival probability of extended SAA(sort whole planning horizon) in " + df.format(sampleNum) + " samples is: " + nf.format(result1[0]));
			error  = 1.96 * Math.sqrt(result1[1]*(1 - result1[1]) / sampleNum);
			thisServiceRate = 1-result1[1];
			System.out.println("final simulated service rate of extended SAA(sort whole planning horizon) " + " is: " + nf.format(thisServiceRate) + " with error " + nf.format(error));
			System.out.println();
			
			double extendSimObj = result1[0];
			double extendSimServiceRate = thisServiceRate;
			double extendSimError = error;
			double extendSimTime = time1;
					
			/*******************************************************************
			 * solve the problem by SDP when there is individual chance constraint approximation,
			 * to get a lower bound
			 */			
			// feasible actions 2
			Function<CashState, double[]> getFeasibleAction2 = s -> {
				int t = s.getPeriod() - 1;
				double thisPeriodServRate = (1 - serviceRate) / T;
				double minQ = Math.ceil(distributions[t].inverseF(1 - thisPeriodServRate)); // minimum ordering quantity in each period 
				return DoubleStream.iterate(minQ, i -> i + stepSize).limit((int) maxOrderQuantity + 1).toArray();
			};
			
			/*******************************************************************
			 * Solve
			 */
			CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction2, stateTransition,
					immediateValue, discountFactor);
			period = 1;		
			initialState = new CashState(period, iniI, iniCash);
			currTime = System.currentTimeMillis();
			recursion.setTreeMapCacheAction();
			double finalValue = recursion.getSurvProb(initialState);
			System.out.println("**********************************************");
			System.out.println("result of SDP with service rate constraint is: ");
			finalValue = finalValue * 100;
			System.out.printf("survival probability for this initial state is: %.2f%%\n", finalValue);
			System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
			time1 = (System.currentTimeMillis() - currTime) / 1000;
			System.out.println("running time is " + time1 + "s");
			double lowerBound = finalValue / 100.0;
			
			
			/*******************************************************************
			 * 
			 * output results to excel
			 */
			System.out.println("");
			double[] out = new double[]{m, serviceRate, sampleNumTotal, iniCash, prices[0], variCostUnits[0], saaObj, saaTime, saaSimObj, saaSimServiceRate, saaSimTime, extendObj, extendTime, extendSimObj, extendSimServiceRate, extendSimTime, staUpperBound, lowerBound, sampleNum, sigmasCoe,holdCostUnit};
			wr.writeToExcelAppend(out, fileName);		 
		}
			    
//			/**
//			 * solve the problem by rolling horizon of extended SAA
//			 * 
//			 */
//			sampleNum = 100; // number of scenarios for rolling SAA
//		    currTime = System.currentTimeMillis();
//		    System.out.println("**********************************************");
//		    System.out.println("after rolling horizon for length " + rollingLength +", result of SAA-scenario tree: ");
//		    double[] resultRolling; 
//		    double[][] scenariosRolling = sampling.generateLHSamples(distributions, sampleNums); //include scenarios in the whole planning horizon
//		    resultRolling = simulation1.rollingHoirzonFurtherExtendSAA(rollingLength, initialState, rollingServiceRate, sampleNumsRolling, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenariosRolling, sampleNum);
//		    time1 = (System.currentTimeMillis() - currTime) / 1000.00;	     
//		    System.out.println("running time is " + time1 + "s");
//		    System.out.println("final simulated survival probability of rolling further SAA in " + df.format(sampleNum) + " samples is: " + nf.format(resultRolling[0]));
//		    double sigma2 = Math.sqrt(resultRolling[1]*(1 - resultRolling[1])/sampleNum);
//			double error2  = 1.96*sigma2;
//			double serviceRate2 = 1 - resultRolling[1];
//			System.out.printf("the service rate for simulated extended SAA rolling horizon is %.4f, with error %.4f\n", serviceRate2, error2);	
		}
	}
}