package sdp.cash;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 14, 2018---10:31:11 AM
*@description:  simulate cash stochastic inventory problem for the algorithm in Chao(2008)'s paper.
*/

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentSkipListMap;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.cash.CashStateXR;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.stat.Tally;


public class CashSimulationXR {

	protected int sampleNum;
	protected Distribution[] distributions;
	CashRecursionXR recursion;
	protected StateTransitionFunction<CashStateXR, Double, Double, CashStateXR> stateTransition; 
	protected ImmediateValueFunction<CashStateXR, Double, Double, Double> immediateValue; 
	double discountFactor;
	double fixOrderCost;
	double price;
	double variOrderCost;
	double holdCost;
	double salvageValue;
	double overheadCost;
	
	Map<CashStateXR, Double> cacheC1Values = new ConcurrentSkipListMap<>();
	Map<CashStateXR, Double> cacheC2Values = new ConcurrentSkipListMap<>();

	/**
	 * simulation for cash sdp problem
	 * @param distributions
	 * @param sampleNum
	 * @param recursion
	 */
	public CashSimulationXR(Distribution[] distributions, int sampleNum,
			CashRecursionXR recursionXR, double discountFactor, double fixOrderCost, double price, 
			double variOrderCost, double holdCost, double salvageValue) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;
		this.recursion = recursionXR;
		this.stateTransition = recursion.getStateTransitionFunction();
		this.immediateValue = recursion.getImmediateValueFunction();		
		this.discountFactor = discountFactor;
		this.fixOrderCost = fixOrderCost;
		this.price = price;
		this.variOrderCost = variOrderCost;
		this.holdCost = holdCost;
		this.salvageValue = salvageValue;
	}
	
	
	public CashSimulationXR(Distribution[] distributions, int sampleNum, ImmediateValueFunction<CashStateXR, Double, Double, Double> immediateValue,
			StateTransitionFunction<CashStateXR, Double, Double, CashStateXR> stateTransition,
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
	public double simulateSDPGivenSamplNum(CashStateXR iniState) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashStateXR state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(Math.max(0, samples[i][t])); // integer samples to test sdp
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniR;
		System.out.println("\nfinal simulated expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}

	public double[] simulateSDPwithErrorConfidence(CashStateXR iniState, double error, double confidence) {
		int minRuns = 1000;   int maxRuns = 1000000;
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();

		Tally costTally = new Tally();
		double[] centerAndRadius = new double[2];
		int sampleNumUse = 0;		
		for(int i = 0; i < minRuns || (centerAndRadius[1]>=centerAndRadius[0]*error && i < maxRuns); i++) {
			double[] realizedDemand = sampling.getNextSample(distributions);
			double sum = 0; CashStateXR state = iniState;
			for (int t = 0; t < realizedDemand.length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(Math.max(0, realizedDemand[t])); // integer samples to test sdp
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
		centerAndRadius[0] += iniState.iniR;
		System.out.println(
				"final simulated expected value in " + confidence*100 + "% confidence level is: " + df1.format(centerAndRadius[0]));
		System.out.println("using " + df2.format(sampleNumUse) + " samples, " +  "confidence interval is [" + "-" + df1.format(centerAndRadius[1]) + ", " + df1.format(centerAndRadius[1]) + "]");
		return centerAndRadius;
	}
	
	/**
	 * @param a* in each period
	 * @return simulate results for the policy proposed by Chao (2008)
	 */
	public double simulateAStar(double[] optY, CashStateXR iniState) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simValues = new double[samples.length];
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; CashStateXR state = iniState;
			for (int t = 0; t < samples[0].length; t++) {
				recursion.getExpectedValue(state);
				double thisY = state.iniR > variOrderCost * optY[t] ? optY[t] : state.iniR / variOrderCost;
				double randomDemand = Math.round(Math.max(0, samples[i][t])); // integer samples to test sdp
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, thisY, randomDemand);
				state = stateTransition.apply(state, thisY, randomDemand);
			}
			simValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simValues).sum()/samples.length + iniState.iniR;
		System.out.println("\nfinal simulated expected value for a* policy with " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}

}
