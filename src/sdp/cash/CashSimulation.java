package sdp.cash;



import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

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
	double fixOrderCost;
	double price;
	double variOrderCost;
	double holdCost;
	double salvageValue;
	double overheadCost;
	
	Map<State, Double> cacheC1Values = new TreeMap<>();
	Map<State, Double> cacheC2Values = new TreeMap<>();

	/**
	 * simulation for cash sdp problem
	 * @param distributions
	 * @param sampleNum
	 * @param recursion
	 */
	public CashSimulation(Distribution[] distributions, int sampleNum,
			CashRecursion recursion, double discountFactor, double fixOrderCost, double price, 
			double variOrderCost, double holdCost, double salvageValue) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;
		this.recursion = recursion;
		this.stateTransition = recursion.getStateTransitionFunction();
		this.immediateValue = recursion.getImmediateValueFunction();		
		this.discountFactor = discountFactor;
		this.fixOrderCost = fixOrderCost;
		this.price = price;
		this.variOrderCost = variOrderCost;
		this.holdCost = holdCost;
		this.salvageValue = salvageValue;
	}
	
	
	public CashSimulation(Distribution[] distributions, int sampleNum, ImmediateValueFunction<CashState, Double, Double, Double> immediateValue,
			StateTransitionFunction<CashState, Double, Double, CashState> stateTransition,
			double discountFactor, double fixOrderCost, double price, 
			double variOrderCost, double holdCost, double salvageValue) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;	
		this.discountFactor = discountFactor;
		this.immediateValue = immediateValue;
		this.stateTransition = stateTransition;
		this.fixOrderCost = fixOrderCost;
		this.price = price;
		this.variOrderCost = variOrderCost;
		this.holdCost = holdCost;
		this.salvageValue = salvageValue;
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
		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
		
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashState state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
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
	 * @param iniState
	 * @param error
	 * @param confidence
	 * @return simulate sdp results with error confidence
	 * @date: Apr 23, 2020, 11:43:23 AM 
	 */
	public double[] simulateSDPwithErrorConfidence(CashState iniState, double error, double confidence) {
		int minRuns = 1000;   int maxRuns = 1000000;
		Sampling.resetStartStream();

		Tally costTally = new Tally();
		double[] centerAndRadius = new double[2];
		int sampleNumUse = 0;		
		for(int i = 0; i < minRuns || (centerAndRadius[1]>=centerAndRadius[0]*error && i < maxRuns); i++) {
			double[] realizedDemand = Sampling.getNextSample(distributions);
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
		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
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
		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
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
							optsCS[t][1] = 0;	
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
		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
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
		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
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

		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
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
	 * compute L(y)
	 * @param  y: order-up-to level y,
	 * @param t : period t + 1
	 * @return value of Ly
	 */
	double Ly(double y, int t) {
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
