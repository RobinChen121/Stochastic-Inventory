package sdp.cash;



import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentSkipListMap;

import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.poi.ss.formula.functions.Count;

import milp.LostSaleChance;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.State;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.stat.Tally;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 14, 2018---10:31:11 AM
*@description:  simulate cash stochastic inventory problem
*/
public class CashSimulation {

	protected int sampleNum;
	protected Distribution[] distributions;
	CashRecursion recursion;
	protected StateTransitionFunction<CashState, Double, Double, CashState> stateTransition; 
	protected ImmediateValueFunction<CashState, Double, Double, Double> immediateValue; 
	double discountFactor;
	
	Map<State, Double> cacheC1Values = new ConcurrentSkipListMap<>();
	Map<State, Double> cacheC2Values = new ConcurrentSkipListMap<>();

	/**
	 * simulation for cash sdp problem
	 * @param distributions
	 * @param sampleNum
	 * @param recursion
	 */
	public CashSimulation(Distribution[] distributions, int sampleNum,
			CashRecursion recursion, double discountFactor) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;
		this.recursion = recursion;
		this.stateTransition = recursion.getStateTransitionFunction();
		this.immediateValue = recursion.getImmediateValueFunction();		
		this.discountFactor = discountFactor;
	}
	
	/**
	 * simulation for SAA
	 * @param distributions
	 * @param sampleNum
	 * @param recursion
	 * @param discountFactor
	 */
	public CashSimulation(Distribution[] distributions, int sampleNum, ImmediateValueFunction<CashState, Double, Double, Double> immediateValue,
			StateTransitionFunction<CashState, Double, Double, CashState> stateTransition) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;
		this.immediateValue = immediateValue;
		this.stateTransition = stateTransition;
	}
	
	
	public CashSimulation(Distribution[] distributions, int sampleNum, ImmediateValueFunction<CashState, Double, Double, Double> immediateValue,
			StateTransitionFunction<CashState, Double, Double, CashState> stateTransition,
			double discountFactor) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;	
		this.discountFactor = discountFactor;
		this.immediateValue = immediateValue;
		this.stateTransition = stateTransition;
	}

	public void setSampleNum(int n) {
		this.sampleNum = n;
	}
	
	
	/**
	 * 
	 * @param iniState
	 * @return simulate sdp in a given number of samples
	 */
	public double simulateSDPGivenSamplNum(CashState iniState) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double mean[] = new double[distributions.length];
		for (int i = 0; i < distributions.length; i++) {
			double sum = 0;
			for (int j = 0; j < sampleNum; j++) {
				sum += samples[j][i];
			}
			mean[i] = sum / sampleNum;
		}
		
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				double thisValue = immediateValue.apply(state, optQ, randomDemand);
				sum += Math.pow(discountFactor, t) * thisValue;
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
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
	
	public double[][] gettTArray(double[][] arr, int t) {
		int T = arr.length;
		int m = arr[0].length;
		double[][] result =  new double[T - t][m];
		for (int i = t; i < T; i++) {
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
	
	public Distribution[] gettTArray(Distribution[] arr, int t) {
		int T = arr.length;
		Distribution[] result =  new Distribution[T - t];
		for (int i = t; i < T; i++)
			result[i-t] = arr[i];
		return result;
	}
	
	
	public double[] simulateSAALostRate(CashState iniState, double iniQ,  double serviceRate, 
			int[] sampleNums, double[] prices, double[] variCostUnits, double[] overheadCosts, double salvageValueUnit, double holdCostUnit,
			double[][] scenarios, int sampleNum) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double[] simuValues = new double[samples.length];
		int T = samples[0].length;
		int lostSaleCount = 0;
		for (int i = 0; i < samples.length; i++) {
			CashState state = iniState;
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
					double nextServiceRate = countLostBefore == true ? 0 : thisServiceRate / distributions[t].cdf(optQ + state.getIniInventory()); // very important
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
					thisServiceRate = nextServiceRate;
				}			
			}
		}

		double simFinalValue = 1 - Arrays.stream(simuValues).sum()/(double)samples.length;
		double lostSaleRate = (double) lostSaleCount / (double) sampleNum;		
		double[] result =  {simFinalValue, lostSaleRate};
		return result;
	}
	

	/**
	 * 
	 * @param iniState
	 * @return simulate sdp in a given number of samples with lost sale rate simulation
	 */
	public double[] simulateSDPGivenSamplNumLostRate(CashState iniState, ImmediateValueFunction<CashState, Double, Double, Double> immediateValue2) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
//		double mean[] = new double[distributions.length];
//		for (int i = 0; i < distributions.length; i++) {
//			double sum = 0;
//			for (int j = 0; j < sampleNum; j++) {
//				sum += samples[j][i];
//			}
//			mean[i] = sum / sampleNum;
//		}
		
		double[] simuValues = new double[samples.length];		
		int lostSaleCount = 0;
		for (int i = 0; i < samples.length; i++) {
			CashState state = iniState;
			boolean countBefore = false;
			boolean countBeforeBankrupt = false;
			for (int t = 0; t < samples[0].length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				if (state.getIniInventory() + optQ < randomDemand - 0.1 && countBefore == false) {
					lostSaleCount ++;
					countBefore = true;
				}
				double thisValue = state.iniCash + immediateValue2.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
				if (thisValue < - 0.1 && countBeforeBankrupt == false) {
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
	 * @param iniState
	 * @return bankrupt probability
	 */
	public double simulateDefaultProb(CashState iniState) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double mean[] = new double[distributions.length];
		for (int i = 0; i < distributions.length; i++) {
			double sum = 0;
			for (int j = 0; j < sampleNum; j++) {
				sum += samples[j][i];
			}
			mean[i] = sum / sampleNum;
		}
		
		double[] simuValues = new double[samples.length];	
		double count = 0;
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			boolean recordMinusCash = false;
			for (int t = 0; t < samples[0].length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				double thisValue = immediateValue.apply(state, optQ, randomDemand);
				sum += Math.pow(discountFactor, t) * thisValue;
				state = stateTransition.apply(state, optQ, randomDemand);
				if (sum + iniState.iniCash < 0 && recordMinusCash == false) {
					count++;
					recordMinusCash = true;
				}
					
			}
			simuValues[i] = sum;				
		}
		return count/(double) sampleNum;
	}

	/**
	 * @param iniState
	 * @param error
	 * @param confidence
	 * @return simulate sdp results with error confidence
	 * @date: Apr 23, 2020, 11:43:23 AM 
	 */
	public double[] simulateSDPwithErrorConfidence(CashState iniState, double error, double confidence) {
		int minRuns = 1000;   int maxRuns = 1000000;
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();		
		Tally costTally = new Tally();
		double[] centerAndRadius = new double[2];
		int sampleNumUse = 0;		
		for(int i = 0; i < minRuns || (centerAndRadius[1]>=centerAndRadius[0]*error && i < maxRuns); i++) {
			double[] realizedDemand = sampling.getNextSample(distributions);
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < realizedDemand.length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(realizedDemand[t]); // integer samples to test sdp
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			costTally.add(sum);
			if(i >= minRuns) 
				costTally.confidenceIntervalNormal(confidence, centerAndRadius);	
			sampleNumUse = i + 1;
		}

		DecimalFormat df1 = new DecimalFormat("0.0000");
		DecimalFormat df2 = new DecimalFormat("###,###");
		centerAndRadius[0] += iniState.iniCash;
		System.out.println(
				"final simulated expected value in " + confidence*100 + "% confidence level is: " + df1.format(centerAndRadius[0]));
		System.out.println("using " + df2.format(sampleNumUse) + " samples, " +  "confidence interval is [" + "-" + df1.format(centerAndRadius[1]) + ", " + df1.format(centerAndRadius[1]) + "]");
		return centerAndRadius;
	}
	
	/**
	 * 
	 * @param iniState
	 * @return simulate (s, C1, C2, S) policy in strong cash constraint
	 */
	public double simulatesCS(CashState iniState, double[][] optsCS, Map<State, Double> cacheC1Values, 
			Map<State, Double> cacheC2Values, double minCashRequired, Double maxQ, double fixOrderCost, double variCost) {
		Sampling.resetStartStream();
		this.cacheC1Values = cacheC1Values;
		this.cacheC2Values = cacheC2Values;
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];		
		int M = 10000;
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t == 0) 
					optQ = iniState.getIniInventory() < optsCS[t][0] ? optsCS[t][3] - iniState.getIniInventory() : 0;
				else {
					double maxOrderQuantity = (int) Math.max(0, (state.iniCash - minCashRequired - fixOrderCost)/variCost);
					maxOrderQuantity = Math.min(maxOrderQuantity, maxQ);
					if (state.getIniInventory() < optsCS[t][0]) {
						if (cacheC1Values.get(new State(state.getPeriod(), state.getIniInventory())) == null)
							optsCS[t][1] = 0;	
						else 
							optsCS[t][1] = cacheC1Values.get(new State(state.getPeriod(), state.getIniInventory()));
						if (cacheC2Values.get(new State(state.getPeriod(), state.getIniInventory())) == null)
							optsCS[t][2] = M;	
						else 
							optsCS[t][2] = cacheC2Values.get(new State(state.getPeriod(), state.getIniInventory()));
					}
					// not include equal
					if (state.getIniInventory() < optsCS[t][0] && state.getIniCash() > optsCS[t][1]
							       && state.getIniCash() < optsCS[t][2])
						optQ = Math.min(maxOrderQuantity, optsCS[t][3] - state.getIniInventory());
					else
						optQ = 0;
				}
				double randomDemand = samples[i][t];
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated (s, C1, C2, S) policy expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}
	
	
	/**
	 * 
	 * @param iniState
	 * @return simulate (s, C1, S) policy in strong cash constraint
	 */
	public double simulatesCS(CashState iniState, double[][] optsCS, Map<State, Double> cacheC1Values, 
			double overheadCost, Double maxQ, double fixOrderCost, double variCost) {
		Sampling.resetStartStream();
		this.cacheC1Values = cacheC1Values;
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];		
		int M = 10000;
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t == 0) {
					optQ = iniState.getIniInventory() < optsCS[t][0] ? optsCS[t][2] - iniState.getIniInventory() : 0;
					double maxOrderQuantity = Math.max(0, (state.iniCash - overheadCost - fixOrderCost)/variCost);
					optQ = Math.min(optQ, maxOrderQuantity);
				}
				else {
					double maxOrderQuantity = Math.max(0, (state.iniCash - overheadCost - fixOrderCost)/variCost);
					maxOrderQuantity = Math.min(maxOrderQuantity, maxQ);
					if (state.getIniInventory() < optsCS[t][0]) {
						if (cacheC1Values.get(new State(state.getPeriod(), state.getIniInventory())) == null)
							optsCS[t][1] = 0;	// 不存在 C(x) 就默认不存资金节点
						else 
							optsCS[t][1] = cacheC1Values.get(new State(state.getPeriod(), state.getIniInventory()));
					}
					// not include equal
					if (state.getIniInventory() < optsCS[t][0] && state.getIniCash() > optsCS[t][1])
						optQ = Math.min(maxOrderQuantity, optsCS[t][2] - state.getIniInventory());
					else
						optQ = 0;
				}
				double randomDemand = samples[i][t];
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated (s, C, S) policy expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}
	
	
	/**
	 * 
	 * @param iniState
	 * @return simulate (s, C, S) policy in strong cash constraint
	 */
	public double simulatesMeanCS(CashState iniState, double[][] optsCS, 
			double minCashRequired, Double maxQ, double fixOrderCost, double variCost) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t == 0) 
					optQ = iniState.getIniInventory() < optsCS[t][0] ? optsCS[t][2] - iniState.getIniInventory() : 0;
				else {
					double maxOrderQuantity = Math.max(0, (state.iniCash - minCashRequired - fixOrderCost)/variCost);
					maxOrderQuantity = Math.min(maxOrderQuantity, maxQ);
					
					// not include equal
					if (state.getIniInventory() < optsCS[t][0] && state.getIniCash() > optsCS[t][1])
						optQ = Math.min(maxOrderQuantity, optsCS[t][2] - state.getIniInventory());
					else
						optQ = 0;
				}
				double randomDemand = samples[i][t];
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated (s, C, S) policy expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}
	
	
	/**
	 * 
	 * @param iniState
	 * @return simulate (s, C, S) policy in overdraft
	 */
	public double simulatesCSDraft(CashState iniState, double[][] optsCS, double minCashRequired, Double maxQ, double fixOrderCost, double variCost) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t == 0) 
					optQ = optsCS[t][2] - optsCS[t][0];
				else {
					double maxOrderQuantity = Math.max(0, (state.iniCash - minCashRequired - fixOrderCost)/variCost);
					maxOrderQuantity = Math.min(maxOrderQuantity, maxQ);
					// not include equal
					if (state.getIniInventory() < optsCS[t][0] && state.getIniCash() > optsCS[t][1])
						optQ = Math.min(maxOrderQuantity, optsCS[t][2] - state.getIniInventory());
					else
						optQ = 0;
				}
				double randomDemand = samples[i][t];
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated (s, C, S) policy expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}
	
	/**
	 * 
	 * @param iniState
	 * @return simulate (s, S) policy in overdraft
	 */
	public double simulatesSOD(CashState iniState, double[][] optsS, double minCashRequired, Double maxQ, double fixOrderCost, double variCost) {
		Sampling.resetStartStream();
		
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t== 0) 
					optQ = optsS[t][1];
				else {
					double maxOrderQuantity = Math.max(0, (state.iniCash -minCashRequired - fixOrderCost)/variCost);
					maxOrderQuantity = Math.min(maxOrderQuantity, maxQ);
					if (state.getIniInventory() < optsS[t][0])
						optQ = Math.min(maxOrderQuantity, optsS[t][1] - state.getIniInventory());
					else
						optQ = 0;
				}
				double randomDemand = samples[i][t];
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated (s, S) policy expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}
	
	/**
	 * 
	 * @param 
	 * @return simulate (s, C, S1, S2) policy in overdraft
	 */
	public double simulatesSOD2(CashState iniState, double[][] optsCS) {
		Sampling.resetStartStream();
		
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t== 0) 
					optQ = optsCS[t][2] - iniState.getIniInventory();
				else {
					if (state.getIniInventory() < optsCS[t][0])
						if (state.getIniCash() <= optsCS[t][1] + 0.1)
							optQ = optsCS[t][2] - state.getIniInventory();
						else 
							optQ = optsCS[t][3] - state.getIniInventory();
					else
						optQ = 0;
				}
				double randomDemand = samples[i][t];
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated (s, C, S1, S2) policy expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}

	/**
	 * compute L(y)
	 * @param  y: order-up-to level y,
	 * @param t : period t + 1
	 * @return value of Ly
	 */
	double Ly(double y, int t, double price, double variOrderCost, double holdCost, double salvageValue) {
		Distribution distribution = distributions[t];
		
		double meanI = 0;
		for (int i = 0; i < y; i++)
			meanI += (y - i) * (distribution.cdf(i + 0.5) - distribution.cdf(i - 0.5));
		
		double Ly = 0;
		if (t == distributions.length - 1)
			Ly = (price - variOrderCost) * y- (price + holdCost - salvageValue) * meanI;
		else
			Ly = (price - variOrderCost) * y- (price + holdCost) * meanI;

		return Ly;
	}
}
