package sdp.cash;

import java.util.Arrays;
import java.util.stream.IntStream;

import milp.LostSaleChance;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;

/**
*@author: zhenchen
*@date: Jun 11, 2023, 12:14:15 PM
*@desp: simulation class for survival maximization problems
*
*/

public class RiskSimulation {
	
	protected int sampleNum;
	protected Distribution[] distributions;
	RiskRecursion recursion;
	protected StateTransitionFunction<RiskState, Double, Double, RiskState> stateTransition; 
	protected ImmediateValueFunction<RiskState, Double, Double, Double> immediateValue;
	
	public RiskSimulation(Distribution[] distributions, int sampleNum, ImmediateValueFunction<RiskState, Double, Double, Double> immediateValue,
			StateTransitionFunction<RiskState, Double, Double, RiskState> stateTransition) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;	
		this.immediateValue = immediateValue;
		this.stateTransition = stateTransition;
	}
	
	public RiskSimulation(Distribution[] distributions, int sampleNum, RiskRecursion recursion) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;
		this.recursion = recursion;
		this.stateTransition = recursion.getStateTransitionFunction();
		this.immediateValue = recursion.getImmediateValueFunction();		
	}
	
	
	public double[] simulateSAA(RiskState iniState, double iniQ,  double serviceRate, 
			int[] sampleNums, double[] prices, double[] variCostUnits, double[] overheadCosts, double salvageValueUnit, double holdCostUnit,
			double[][] scenarios, int sampleNum) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double[] simuValues = new double[samples.length];
		int T = samples[0].length;
		int lostSaleCount = 0;
		for (int i = 0; i < samples.length; i++) {
			RiskState state = iniState;
			boolean countBeforeBankrupt = false;
			boolean countLostBefore = false;
			double thisValue = 0;
			double optQ = 0;
			double thisServiceRate = 1;
			for (int t = 0; t < T; t++) {
				if (t == 0) {
					optQ = iniQ;
					thisServiceRate = serviceRate;
				}				
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				thisValue = state.iniCash + immediateValue.apply(state, optQ, randomDemand);
				if (state.getIniInventory() + optQ < randomDemand - 0.1 && countLostBefore == false) {
					lostSaleCount ++;
					countLostBefore = true;
				}
				if (thisValue < - 0.1 && countBeforeBankrupt == false) {
					simuValues[i] = 1;
					countBeforeBankrupt = true;
				}
				if (t < T - 1) {										
					double nextServiceRate = 0;
					
//					double thisPeriodServRate = distributions[t].cdf(optQ + state.getIniInventory());
//					nextServiceRate = thisPeriodServRate < thisServiceRate ? thisServiceRate : thisServiceRate / thisPeriodServRate;				
					
					double meanDemandSum = IntStream.range(0, T).mapToDouble(j -> distributions[j].getMean()).sum();
					double rollingDemandSum = IntStream.range(t, T).mapToDouble(j -> distributions[j].getMean()).sum();
					double portion = rollingDemandSum / meanDemandSum;
					nextServiceRate = Math.pow(serviceRate, portion);
					
					state = stateTransition.apply(state, optQ, randomDemand);
					double iniCash = state.iniCash;
					double iniI = state.getIniInventory();
					double[] nextPrices = gettTArray(prices, t + 1);
					int[] nextSampleNums = gettTArray(sampleNums, t + 1);
					double[] nextVariCostUnits = gettTArray(variCostUnits, t + 1);
					Distribution[] nexDistributions = gettTArray(distributions, t + 1);
					double[] nextOverheadCosts = gettTArray(overheadCosts, t + 1);
					double[][] nextScenarios = gettTArray(scenarios, t + 1);
					LostSaleChance model = new LostSaleChance(nexDistributions, nextSampleNums, iniCash, iniI, nextPrices, nextVariCostUnits, 
							salvageValueUnit, holdCostUnit, nextOverheadCosts, nextServiceRate, nextScenarios);
					double[] result = model.solveMaxSurvival();	
					optQ = result[0];
					if (state.iniCash < 0)
						optQ = 0;
					thisServiceRate = nextServiceRate;
				}			
			}
		}

		double simObj = 1 - Arrays.stream(simuValues).sum()/(double)samples.length;
		double simLostRate = (double) lostSaleCount / (double) sampleNum;		
		int n = samples.length;
		double sigma = Math.sqrt(simObj*(1 - simObj)/n);
		double lowCI = simObj - 1.96*sigma;
		double upCI = simObj + 1.96*sigma;
		double error = 1.96*sigma;
		System.out.printf("the confidence interval for simulated  SAA objective is [%.4f, %.4f], with error is %.4f. \n", lowCI, upCI, error);	
		double[] result =  {simObj, simLostRate};
		return result;
	}
	
	public double[] simulateScenarioTree(RiskState iniState, double iniQ,  double serviceRate, 
			int[] sampleNums, double[] prices, double[] variCostUnits, double[] overheadCosts, double salvageValueUnit, double holdCostUnit,
			double[][] scenarios, int sampleNum) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double[] simuValues = new double[samples.length];
		int T = samples[0].length;
		int lostSaleCount = 0;
		for (int i = 0; i < samples.length; i++) {
			RiskState state = iniState;
			boolean countBeforeBankrupt = false;
			boolean countLostBefore = false;
			double thisValue = 0;
			double optQ = 0;
			double thisServiceRate = 1;
			for (int t = 0; t < T; t++) {
				if (t == 0) {
					optQ = iniQ;
					thisServiceRate = serviceRate;
				}				
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				thisValue = state.iniCash + immediateValue.apply(state, optQ, randomDemand);
				if (state.getIniInventory() + optQ < randomDemand - 0.1 && countLostBefore == false) {
					lostSaleCount ++;
					countLostBefore = true;
				}
				if (thisValue < - 0.1 && countBeforeBankrupt == false) {
					simuValues[i] = 1;
					countBeforeBankrupt = true;
				}
				if (t < T - 1) {					
					double nextServiceRate = 0;
					
//                  double thisPeriodServRate = distributions[t].cdf(optQ + state.getIniInventory());	
//					nextServiceRate = thisPeriodServRate < thisServiceRate ? thisServiceRate : thisServiceRate / thisPeriodServRate;

	                double meanDemandSum = IntStream.range(0, T).mapToDouble(j -> distributions[j].getMean()).sum();
					double rollingDemandSum = IntStream.range(t, T).mapToDouble(j -> distributions[j].getMean()).sum();
					double portion = rollingDemandSum / meanDemandSum;
					nextServiceRate = Math.pow(serviceRate, portion);
//					double nextServiceRate2 = Math.pow(serviceRate, portion);
//					nextServiceRate = Math.max(nextServiceRate, nextServiceRate2);
					
//					if (countLostBefore == true)
//						nextServiceRate = 0;
//					else {
//						nextServiceRate = thisPeriodServRate < thisServiceRate ? thisServiceRate : thisServiceRate / thisPeriodServRate;
//					} // very important
					state = stateTransition.apply(state, optQ, randomDemand);
					double iniCash = state.iniCash;
					double iniI = state.getIniInventory();
					double[] nextPrices = gettTArray(prices, t + 1);
					int[] nextSampleNums = gettTArray(sampleNums, t + 1);
					double[] nextVariCostUnits = gettTArray(variCostUnits, t + 1);
					Distribution[] nexDistributions = gettTArray(distributions, t + 1);
					double[] nextOverheadCosts = gettTArray(overheadCosts, t + 1);
					double[][] nextScenarios = gettTArray(scenarios, t + 1);
					LostSaleChance model = new LostSaleChance(nexDistributions, nextSampleNums, iniCash, iniI, nextPrices, nextVariCostUnits, 
							salvageValueUnit, holdCostUnit, nextOverheadCosts, nextServiceRate, nextScenarios);
					double[] result = model.solveScenario();	
					optQ = result[0];
					if (state.iniCash < 0)
						optQ = 0;
					thisServiceRate = nextServiceRate;
				}			
			}
		}

		int n = samples.length;
		double simFinalValue = 1 - Arrays.stream(simuValues).sum()/(double)samples.length;
		double sigma = Math.sqrt(simFinalValue*(1 - simFinalValue)/n);
		double lowCI = simFinalValue - 1.96*sigma;
		double upCI = simFinalValue + 1.96*sigma;
		double error  = 1.96*sigma;
		System.out.printf("the confidence interval for simulated scenario tree is [%.4f, %.4f], with error %.4f\n", lowCI, upCI, error);
		double lostSaleRate = (double) lostSaleCount / (double) sampleNum;		
		double[] result =  {simFinalValue, lostSaleRate};
		return result;
	}
	
	/**
	 * 
	 * @param iniState
	 * @return simulate sdp in a given number of samples with lost sale rate simulation
	 */
	public double[] simulateLostSale(RiskState iniState, ImmediateValueFunction<RiskState, Double, Double, Double> immediateValue2) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double[] simuValues = new double[samples.length];		
		int lostSaleCount = 0;
		for (int i = 0; i < samples.length; i++) {
			RiskState state = iniState;
			boolean countBefore = false;
			boolean countBeforeBankrupt = state.getBankruptBefore();
			for (int t = 0; t < samples[0].length; t++)
			{
				recursion.getSurvProb(state);
				double optQ = recursion.getAction(state);
				if (state.iniCash < 0)
					optQ = 0;
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				if (state.getIniInventory() + optQ < randomDemand && countBefore == false) {
					lostSaleCount ++;
					countBefore = true;
				}
				double thisValue = state.iniCash + immediateValue2.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
				if (thisValue < 0 && countBeforeBankrupt == false) {					
					simuValues[i] = 1;
					countBeforeBankrupt = true;
				}
			}
			
		}
		double simFinalValue = 1 - Arrays.stream(simuValues).sum()/(double)samples.length;
		double lostSaleRate = (double) lostSaleCount / (double) sampleNum;		
		double[] result =  {simFinalValue, lostSaleRate};
		return result;
	}
	
	/**
	 * rolling horizon for the further extended SAA
	 * @param r is the rolling length
	 * @param iniState
	 * @param iniQ
	 * @param serviceRate
	 * @param sampleNums
	 * @param prices
	 * @param variCostUnits
	 * @param overheadCosts
	 * @param salvageValueUnit
	 * @param holdCostUnit
	 * @param scenarios
	 * @param sampleNum
	 * @return
	 */
	public double[] rollingHoirzon(int r, RiskState iniState, double serviceRate, 
			int[] sampleNums, double[] prices, double[] variCostUnits, double[] overheadCosts, double salvageValueUnit, double holdCostUnit,
			double[][] scenarios, int sampleNum) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double iniCash = iniState.iniCash;
		double iniI = iniState.getIniInventory();
		double[] nextPrices = gettTArray(prices, 0,r);
		int[] nextSampleNums = gettTArray(sampleNums, 0, r);
		double[] nextVariCostUnits = gettTArray(variCostUnits, 0, r);
		Distribution[] nexDistributions = gettTArray(distributions, 0, r);
		double[] nextOverheadCosts = gettTArray(overheadCosts, 0, r);
		double[][] nextScenarios = gettTArray(scenarios, 0, r);
		LostSaleChance model = new LostSaleChance(nexDistributions, nextSampleNums, iniCash, iniI, nextPrices, nextVariCostUnits, 
				salvageValueUnit, holdCostUnit, nextOverheadCosts, serviceRate, nextScenarios);
		
//		double[] result = model.solveSortWhole();
		double[] result = model.solveMaxSurvival();	
//		double[] result = model.solveScenario();	
		double iniQ = result[0];
		
		double[] simuValues = new double[samples.length];
		int T = samples[0].length;
		int lostSaleCount = 0;
		for (int i = 0; i < samples.length; i++) {
			RiskState state = iniState;
			boolean countBeforeBankrupt = false;
			boolean countLostBefore = false;
			double thisValue = 0;
			double optQ = 0;
			double thisServiceRate = 1;
			for (int t = 0; t < T; t++) {
				if (t == 0) {
					optQ = iniQ;
					thisServiceRate = serviceRate;
				}				
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				thisValue = state.iniCash + immediateValue.apply(state, optQ, randomDemand);
				if (state.getIniInventory() + optQ < randomDemand - 0.1 && countLostBefore == false) {
					lostSaleCount ++;
					countLostBefore = true;
				}
				if (thisValue < - 0.1 && countBeforeBankrupt == false) {
					simuValues[i] = 1;
					countBeforeBankrupt = true;
				}
				if (t < T - 1) {									
					double nextServiceRate;
										
					double meanDemandSum = IntStream.range(0, T).mapToDouble(j -> distributions[j].getMean()).sum();
					double rollingDemandSum = IntStream.range(t, Math.min(t+r, T)).mapToDouble(j -> distributions[j].getMean()).sum();
					double portion = rollingDemandSum / meanDemandSum;
					nextServiceRate = Math.pow(serviceRate, portion);
					
//					double thisPeriodServRate = distributions[t].cdf(optQ + state.getIniInventory());
//					nextServiceRate = thisPeriodServRate > nextServiceRate ? thisPeriodServRate : nextServiceRate;
					
					state = stateTransition.apply(state, optQ, randomDemand);
					iniCash = state.iniCash;
					iniI = state.getIniInventory();
					nextPrices = gettTArray(prices, t + 1, r);
					nextSampleNums = gettTArray(sampleNums, t + 1, r);

					nextVariCostUnits = gettTArray(variCostUnits, t + 1, r);
					nexDistributions = gettTArray(distributions, t + 1, r);
					nextOverheadCosts = gettTArray(overheadCosts, t + 1, r);
					nextScenarios = gettTArray(scenarios, t + 1, r);
					model = new LostSaleChance(nexDistributions, nextSampleNums, iniCash, iniI, nextPrices, nextVariCostUnits, 
							salvageValueUnit, holdCostUnit, nextOverheadCosts, nextServiceRate, nextScenarios);

//					result = model.solveSortWhole();
					result = model.solveMaxSurvival();
//					result = model.solveScenario();
					optQ = result[0];
					if (countBeforeBankrupt == true)
						optQ = 0;
					thisServiceRate = nextServiceRate;
				}			
			}
		}

		int n = samples.length;
		double simFinalValue = 1 - Arrays.stream(simuValues).sum()/(double)samples.length;
		double sigma = Math.sqrt(simFinalValue*(1 - simFinalValue)/n);
		double lowCI = simFinalValue - 1.96*sigma;
		double upCI = simFinalValue + 1.96*sigma;
		double error  = 1.96*sigma;
		System.out.printf("the confidence interval for simulated rolling horizon is [%.4f, %.4f], with error %.4f\n", lowCI, upCI, error);
		double lostSaleRate = (double) lostSaleCount / (double) sampleNum;	
		double[] results =  {simFinalValue, lostSaleRate, iniQ};
		return results;
	}
		
	
	/**
	 * @param arr
	 * @param t
	 * @return an array from index t to T-1
	 */
	public double[] gettTArray(double[] arr, int t) {
		int T = arr.length;
		double[] result =  new double[T - t];
		for (int i = t; i < T; i++)
			result[i-t] = arr[i];
		return result;
	}
	
	/**
	 * @param arr
	 * @param t
	 * @return an array from index t to R (R = t + r > T ? T: t + r; )
	 */
	public double[] gettTArray(double[] arr, int t, int r) {
		int T = arr.length;
		int R = t + r > T ? T: t + r; 
		double[] result =  new double[R - t];
		for (int i = t; i < R; i++)
			result[i-t] = arr[i];
		return result;
	}
	
	public double[][] gettTArray(double[][] arr, int t) {
		int T = arr.length;
		
		double[][] result =  new double[T - t][];
		for (int i = t; i < T; i++) {
			int m = arr[i].length;
			result[i-t] = new double[m];
			for (int j = 0; j < m; j++) 		
				result[i-t][j] = arr[i][j];		
		}
		return result;
	}
	
	public double[][] gettTArray2(double[][] arr, int t) {
		int T = arr[0].length;
		int m = arr.length;
		
		double[][] result =  new double[m][T-t];
		for (int i = t; i < T; i++) {
			for (int j = 0; j < m; j++) 		
				result[j][i-t] = arr[j][i];		
		}
		return result;
	}
	
	public double[][] gettTArray(double[][] arr, int t, int r) {
		int T = arr.length;
		int R = t + r > T ? T: t + r; 
		double[][] result =  new double[R - t][];
		for (int i = t; i < R; i++) {
			int m = arr[i].length;
			result[i-t] = new double[m];
			for (int j = 0; j < m; j++) 		
				result[i-t][j] = arr[i][j];		
		}
		return result;
	}
	
	public int[] gettTArray(int[] arr, int t) {
		int T = arr.length;
		int[] result =  new int[T - t];
		for (int i = t; i < T; i++)
			result[i-t] = arr[i];
		return result;
	}
	
	public int[] gettTArray(int[] arr, int t, int r) {
		int T = arr.length;
		int R = t + r > T ? T: t + r; 
		int[] result =  new int[R - t];
		for (int i = t; i < R; i++)
			result[i-t] = arr[i];
		return result;
	}
	
	public Distribution[] gettTArray(Distribution[] arr, int t) {
		int T = arr.length;
		Distribution[] result =  new Distribution[T - t];
		for (int i = t; i < T; i++)
			result[i-t] = arr[i];
		return result;
	}
	
	public Distribution[] gettTArray(Distribution[] arr, int t, int r) {
		int T = arr.length;
		int R = t + r > T ? T: t + r; 
		Distribution[] result =  new Distribution[R - t];
		for (int i = t; i < R; i++)
			result[i-t] = arr[i];
		return result;
	}

}


