package sdp.chance;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import milp.LostSaleChance;
import sdp.cash.CashRecursion;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;

/*
* @author chen
* @date 2022 Nov 2, 11:27:52
* @describe: test the simulated lost sale rate for Tom's model (2002)
*
*/

public class ChanceTestingLostsaleRate {

	public static void main(String[] args) {
		WriteToExcelTxt wr = new WriteToExcelTxt();
		String fileName = "JointChanceSAA.xls";
		String headString =  
				"demand mode" + "\t" + "rolling length" + "\t"  + "serviceRate" + "\t" + "sample one period" +  "\t"+ "scenario number" + "\t" + "iniCash" + "\t" + "price" + "\t" + "variCost" + "\t" + "rolling time" + "\t" + "rolling obj" + "\t" +
		         "rolling service rate" + "\t" + "sigmasCoe" + "\t" +"holdingCost"+ "\t" + "Tom obj"
					+ "\t" + "sim Tom service" + "\t" + "individual chance obj" + "\t" + "sim individual chance service";
		wr.writeToFile(fileName, headString);
		
		double[][] meanDemands = {{30,30,30,30},
				{49,50,46,38},
				{14, 23, 33, 42},
				{30,6,30,54},
				{30,21,30,39},
				{27,24,23,35},
				{15,140,147,74},
				{24,118,86,117},
				{35,43,59,55},
				{56,84,67,155}
		};
		
		double iniI = 0;
		double trunQuantile = 0.99;
		
		int T = meanDemands[0].length;
		int[] sampleNums = new int[T];
		double[] prices = new double[T];
		double[] variCostUnits = new double[T];
		double[] overheadCosts = new double[T];	
		double sigmasCoe = 0.25; // for use in normal distribution
		
		Arrays.fill(prices, 5);	//
		Arrays.fill(variCostUnits, 1); //
		Arrays.fill(overheadCosts, 80); // overhead costs
		
		double holdCostUnit = 0.5;
		double salvageValueUnit = 0.5;
		double maxOrderQuantity = 300; // maximum ordering quantity when having enough cash
		
		int[] sampleNumsRolling = new int[T];
		int sampleOnePeriod = 9; 
		Arrays.fill(sampleNumsRolling,sampleOnePeriod);
		
		double serviceRate = 0.6; // the higher value results in slower running speed. maximum negative possibility rate is 1 - serviceRate. 
		int etaFrac = 2;
		double serviceRateSAA = 1 - (1- serviceRate)/etaFrac;
		int rollingLength = 3; // rolling horizon length	
		
		double iniCash = 60;
		for (int m = 7; m < 10; m++) {
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
						
				int sampleNum = 100; // number of scenarios for rolling SAA
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
				
				
				double finalValue1 = 0;
				double[] result1 = new double[2];
				double finalValue2 = 0;
				double[] result2 = new double[2];
				
				if (kk < 1) {
				/**
				 * solve the problem by SDP when there is no joint chance constraint
				 */	
				// feasible actions
			    Function<CashState, double[]> getFeasibleAction1 = s -> {
					return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxOrderQuantity + 1).toArray();
				};
				
				
				/*******************************************************************
				 * Solve
				 */
				CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction1, stateTransition,
						immediateValue, discountFactor);
				
				double time;
				
//				currTime = System.currentTimeMillis();
//				recursion.setTreeMapCacheAction();
//				finalValue1 = recursion.getSurvProb(initialState);
//				System.out.println("**********************************************");
//				System.out.println("result of SDP with no service rate constraint is: ");
//				System.out.println("survival probability for this initial state is: " + nf.format(finalValue1));
//				System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
//				time = (System.currentTimeMillis() - currTime) / 1000;
//				System.out.println("running time is " + time + "s");
				
				/*******************************************************************
				 * Simulating sdp results
				 */
				int sampleNumSim = 1000;
				CashSimulation simulation = new CashSimulation(distributions, sampleNumSim, recursion, discountFactor); // no need to add overheadCost in this class
				
//				result1 = simulation.simulateSDPGivenSamplNum(initialState, immediateValue);
//				System.out.println("final simulated survival probability in " + df.format(sampleNumSim) + " samples is: " + nf.format(result1[0]));
//				System.out.println("final simulated service rate " + " is: " + nf.format(1-result1[1])); 
				
				/*******************************************************************
				 * solve the problem by SDP when there is individual chance constraint approximation
				 */			
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
				finalValue2 = recursion.getSurvProb(initialState);
				System.out.println("**********************************************");
				System.out.println("result of SDP with service rate constraint is: ");
				finalValue2 = finalValue2 * 100;
				System.out.printf("survival probability for this initial state is: %.2f%%\n", finalValue2);
				System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
				time = (System.currentTimeMillis() - currTime) / 1000;
				System.out.println("running time is " + time + "s");
				
				/*******************************************************************
				 * Simulating sdp results
				 */		
				sampleNumSim = 1000;
				simulation = new CashSimulation(distributions, sampleNumSim, recursion, discountFactor); // no need to add overheadCost in this class
				result2 = simulation.simulateLostSale(initialState, immediateValue);
				System.out.println("final simulated survival probability in " + df.format(sampleNumSim) + " samples is: " + nf.format(result2[0]));
				System.out.println("final simulated service rate " + " is: " + nf.format(1-result2[1]));
				}
				
				/*******************************************************************
				 * 
				 * output results to excel
				 */
				System.out.println("");
				double[] out = new double[]{m, rollingLength, serviceRate, sampleOnePeriod, sampleNum, iniCash, prices[0], variCostUnits[0], time1, resultRolling[0], serviceRate2, sigmasCoe, holdCostUnit,
												finalValue1, 1-result1[1], finalValue2, 1-result2[1]};
				wr.writeToExcelAppend(out, fileName);
			}
				
		}
	
	}

}
