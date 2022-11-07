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
 * need seasonal price and vari cost for the rolling horizon testing
 *
 */
public class ChanceCashTestingRollingHorizon {
	
	public static void main(String[] args) {
		WriteToExcel wr = new WriteToExcel();
		String fileName = "RollingHorizonPerformance.xls";
		String headString =  
				"demand mode" + "\t" + "rolling length" + "\t"  + "serviceRate" + "\t" + "sample one period" +  "\t"+ "scenario number" + "\t" + "iniCash" + "\t" + "price" + "\t" + "variCost" + "\t" + "rolling time" + "\t" + "rolling obj" + "\t" +
		         "rolling service rate" + "\t" + "sigmasCoe" + "\t" +"holdingCost";
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
						
		double iniI = 0;
		double trunQuantile = 0.99;
		double serviceRate = 0.8; // the higher value results in slower running speed. maximum negative possibility rate is 1 - serviceRate. 
		
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
		Arrays.fill(prices, 5);	//
		Arrays.fill(variCostUnits, 1); //
		Arrays.fill(overheadCosts, 80); // overhead costs
		
		double holdCostUnit = 0.5;
		double salvageValueUnit = 0.5;
		double maxOrderQuantity = 300; // maximum ordering quantity when having enough cash
		
		int[] sampleNumsRolling = new int[T];
		int sampleOnePeriod = 13; 
		Arrays.fill(sampleNumsRolling,sampleOnePeriod);
		
		int rollingLength = 3; // rolling horizon length		
		int instanceNum = meanDemands.length;
		double iniCash = 100;
		for (int m = 2; m < 3; m++) {
			for(int kk = 0; kk < 10; kk++) {	
			/**
			 * solve the problem by rolling horizon of SAA
			 * 
			 */
			
			int instanceIndex = m;
			Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
						.mapToObj(i -> new NormalDist(meanDemands[instanceIndex][i], sigmasCoe*meanDemands[instanceIndex][i])).toArray(Distribution[]::new); //can be changed to other distributions
						//.mapToObj(i -> new PoissonDist(meanDemands[instanceIndex][i])).toArray(Distribution[]::new); // Poisson demand
			
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
					
			int sampleNum = 10; // number of scenarios for rolling SAA
			int period = 1;		
		    CashState initialState = new CashState(period, iniI, iniCash);			    			    			    	    
		    CashSimulation simulation1 = new CashSimulation(distributions, sampleNum, immediateValue, stateTransition); // no need to add overheadCost in this class
		    		
		    long currTime = System.currentTimeMillis();
		    System.out.println("**********************************************");
		    System.out.println("after rolling horizon for length " + rollingLength +", result of SAA-scenario tree: ");
		    double[] resultRolling; 
		    Sampling sampling = new Sampling();
		    
		    double meanDemandSum = Arrays.stream(meanDemands[m]).sum();
			double rollingDemandSum = Arrays.stream(meanDemands[m]).limit(rollingLength).sum();
			double portion = rollingDemandSum / meanDemandSum;
			double rollingServiceRate = Math.pow(serviceRate, portion);
			DecimalFormat df = new DecimalFormat("###, ###");
			NumberFormat nf = NumberFormat.getPercentInstance();
			nf.setMinimumFractionDigits(5);	
			
		    double[][] scenariosRolling = sampling.generateLHSamples(distributions, sampleNumsRolling); //include scenarios in the whole planning horizon
		    resultRolling = simulation1.rollingHoirzonFurtherExtendSAA(rollingLength, initialState, rollingServiceRate, sampleNumsRolling, prices, variCostUnits, overheadCosts, salvageValueUnit, holdCostUnit, scenariosRolling, sampleNum);
		    double time1 = (System.currentTimeMillis() - currTime) / 1000.00;	     
		    System.out.println("running time is " + time1 + "s");
		    System.out.println("final simulated survival probability of rolling further SAA in " + df.format(sampleNum) + "samples is: " + nf.format(resultRolling[0]));
		    double sigma2 = Math.sqrt(resultRolling[1]*(1 - resultRolling[1])/sampleNum);
			double error2  = 1.96*sigma2;
			double serviceRate2 = 1 - resultRolling[1];
			System.out.printf("the service rate for simulated extended SAA rolling horizon is %.4f, with error %.4f\n", serviceRate2, error2);
			
			/*******************************************************************
			 * 
			 * output results to excel
			 */
			System.out.println("");
			double[] out = new double[]{m, rollingLength, serviceRate, sampleOnePeriod, sampleNum, iniCash, prices[0], variCostUnits[0], time1, resultRolling[0], serviceRate2, sigmasCoe, holdCostUnit};
			wr.writeToExcelAppend(out, fileName);
			
		}			
		}
	}
}