package sdp.chance;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

import milp.LostSaleChance;
import milp.LostSaleChanceTesting;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import sdp.write.WriteToExcel;
import umontreal.ssj.probdist.BinomialDist;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;

/**
 * @author chen
 * @email: 15011074486@163.com
 * @Date: 2022.04.19, 10:25:02
 * @Description: test joint chance of sampling with replacement
 * 
 */
public class ChanceCashTestingSampleReplace {
	public static void main(String[] args) {
		WriteToExcel wr = new WriteToExcel();
		String fileName = "JointChanceSAA5Periods.xls";
		String headString =  
				"demand mode" + "\t" + "serviceRate" + "\t" + "scenario number" + "\t" + "iniCash" + "\t" + "price" + "\t" + "variCost" + "\t" + "SAA obj" + "\t" + "time" + "\t" + "sim SAA obj" + "\t" +
		         "sim SAA service rate" + "\t" + "sim saa time" + "\t" + "extend SAA obj" + "\t" + "time" + "\t" + "sim extend SAA obj" + "\t" 
					 + "sim extend SAA service" + "\t" + "sim extend time" + "\t" + "upper bound" + "\t" + "lower bound" +  "\t" + "simSampleNum"
					 + "\t" + "sigmasCoe";
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
		int sampleTotalNumber = 3000;
		int sampleNumSimulate = 200;  // sample number in simulation
		
		int sampleNumPeriod = 10; // number of samples in each period
		Arrays.fill(sampleNums, sampleNumPeriod);
		for (int t = 0; t < T; t++) {
			int m = t % 4;
			prices[t] = seasonalPrice[m];
			variCostUnits[t] = seasonalVariCost[m];
		}
		Arrays.fill(prices, 16);	
		Arrays.fill(variCostUnits, 10);
		Arrays.fill(overheadCosts, 100); // overhead costs
		double holdCostUnit = 1;
		double salvageValueUnit = 5;
		double maxOrderQuantity = 300; // maximum ordering quantity when having enough cash
					
		int instanceNum = 1; //meanDemands.length;
		int runNum = 10;
		Random rd = new Random();
		
		for (int m = 0; m < 1; m++) {
			long currTime = System.currentTimeMillis();
			double time1;
			LostSaleChanceTesting model;
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
			Arrays.fill(sampleNums, sampleNumPeriod);	// each period 10 samples
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			
			int instanceIndex = m;
			Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
						.mapToObj(i -> new NormalDist(meanDemands[instanceIndex][i], sigmasCoe*meanDemands[instanceIndex][i])).toArray(Distribution[]::new); //can be changed to other distributions
						//.mapToObj(i -> new PoissonDist(meanDemands[instanceIndex][i])).toArray(Distribution[]::new); // Poisson demand
		
			Sampling sampling = new Sampling();
			double[][] scenarios = new double[T][];
			for (int k = 0; k < 1; k++) { // run number
				// generate scenarios
				scenarios = sampling.generateLHSamples(distributions, sampleNums);		
				// sample without replacement
				
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
				model = new LostSaleChanceTesting(distributions, demandSamples, iniCash, iniI, prices, variCostUnits, 
									salvageValueUnit, holdCostUnit, overheadCosts, serviceRate);		
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
		    CashSimulation simulation1 = new CashSimulation(distributions, sampleNumSimulate, immediateValue, stateTransition); // no need to add overheadCost in this class
		    double error;
		    double thisServiceRate;		
		    double[] result1;     
		    
			/**
		     * Simulate the restult of SAA
			 */			    			    
		    double[][] demandSamples = new double[sampleTotalNumber][T];
			for (int i = 0; i < sampleTotalNumber; i++) {
				for (int t = 0; t < T; t++) {
					int index = rd.nextInt(sampleNumPeriod);
					demandSamples[i][t] = scenarios[t][index];
				}
			}
		    currTime = System.currentTimeMillis();
		    result1 = simulation1.simulateSAATesting(initialState, resultSAA[0], serviceRate, scenarios, sampleNums, demandSamples, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, sampleNumSimulate);
		    time1 = (System.currentTimeMillis() - currTime) / 1000.00;	
		    System.out.println("running time is " + time1 + "s");
		    System.out.println("final simulated survival probability of SAA in " + df.format(sampleNumSimulate) + " samples is: " + nf.format(result1[0]));
			error  = 1.96 * Math.sqrt(result1[1]*(1 - result1[1]) / sampleNumSimulate);
			thisServiceRate = 1-result1[1];
			System.out.println("final simulated service sale rate of SAA " + " is: " + nf.format(thisServiceRate) + " with error " + nf.format(error)); 
			System.out.println();
			
			double saaSimObj = result1[0];
			double saaSimError = error;
			double saaSimTime = time1;
			double saaSimServiceRate = thisServiceRate;
			
			
						  
		
			
		}			
	}
}
