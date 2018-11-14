package sdp.inventory;

import java.text.DecimalFormat;
import java.util.Arrays;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.stat.Tally;

/**
 *@author: Zhen Chen
 *@email: 15011074486@163.com
 *@date: Jul 11, 2018---4:21:18 PM
 *@description:  Simulation the sdp
 */

public class Simulation {

	protected int sampleNum;
	protected Distribution[] distributions;
	Recursion recursion;
	protected StateTransitionFunction<State, Double, Double, State> stateTransition; 
	protected ImmediateValueFunction<State, Double, Double, Double> immediateValue; 

	/**
	 * simulation for inventory only problem
	 * @param distributions
	 * @param sampleNum
	 * @param recursion
	 */
	public Simulation(Distribution[] distributions, int sampleNum,
			Recursion recursion) {
		this.distributions = distributions;
		this.sampleNum = sampleNum;
		this.recursion = recursion;
		this.stateTransition = recursion.getStateTransitionFunction();
		this.immediateValue = recursion.getImmediateValueFunction();		
	}

	public void setSampleNum(int n) {
		this.sampleNum = n;
	}
	
	
	/**
	 * 
	 * @param iniState
	 * @return simulate sdp in a given number of samples
	 */
	public double simulateSDPGivenSamplNum(State iniState) {
		Sampling.resetStartStream();

		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; State state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
				sum += immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length;
		System.out.println("\nfinal simulated expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}

	public double[] simulateSDPwithErrorConfidence(State iniState, double error, double confidence) {
		int minRuns = 1000;   int maxRuns = 1000000;
		Sampling.resetStartStream();

		Tally costTally = new Tally();
		double[] centerAndRadius = new double[2];
		int sampleNumUse = 0;		
		for(int i = 0; i < minRuns || (centerAndRadius[1]>=centerAndRadius[0]*error && i < maxRuns); i++) {
			double[] realizedDemand = Sampling.getNextSample(distributions);
			double sum = 0; State state = iniState;
			for (int t = 0; t < realizedDemand.length; t++)
			{
				recursion.getExpectedValue(state);
				double optQ = recursion.getAction(state);
				double randomDemand = Math.round(realizedDemand[t]); // integer samples to test sdp
				sum += immediateValue.apply(state, optQ, randomDemand);
				state = stateTransition.apply(state, optQ, randomDemand);
			}
			costTally.add(sum);
			if(i >= minRuns) 
				costTally.confidenceIntervalNormal(confidence, centerAndRadius);	
			sampleNumUse = i + 1;
		}

		DecimalFormat df1 = new DecimalFormat("0.0000");
		DecimalFormat df2 = new DecimalFormat("###,###");
		System.out.println(
				"final simulated expected value in " + confidence*100 + "% confidence level is: " + df1.format(centerAndRadius[0]));
		System.out.println("using " + df2.format(sampleNumUse) + " samples, " +  "confidence interval is [" + "-" + df1.format(centerAndRadius[1]) + ", " + df1.format(centerAndRadius[1]) + "]");
		return centerAndRadius;
	}
	
	
	
}
