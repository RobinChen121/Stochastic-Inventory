/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 25, 2019, 7:15:40 PM
 * @Desc: simulate the sdp results for multi item cash constrained problem
 *
 *
 * 
 */
package sdp.cash.multiItem;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdistmulti.BiNormalDist;

public class CashSimulationMulti {
	int sampleNum;
	BiNormalDist[] distributions;
	
	Distribution[][] distributionsMulti;
	
	double discountFactor;
	
	CashRecursionMulti recursion;
	
	StateTransitionFunction<CashStateMulti, Actions, Demands, CashStateMulti> stateTransition;
	ImmediateValueFunction<CashStateMulti, Actions, Demands, Double> immediateValue;
	
	public CashSimulationMulti(int sampleNum, BiNormalDist[] distributions, double discountFactor, CashRecursionMulti recursion, StateTransitionFunction<CashStateMulti, Actions, Demands, CashStateMulti> stateTransition,
							ImmediateValueFunction<CashStateMulti, Actions, Demands, Double> immediateValue) {
		this.sampleNum = sampleNum;
		this.distributions = distributions;
		this.discountFactor = discountFactor;
		this.recursion = recursion;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
	}
	
	public CashSimulationMulti(int sampleNum, Distribution[][] distributions, double discountFactor, CashRecursionMulti recursion, StateTransitionFunction<CashStateMulti, Actions, Demands, CashStateMulti> stateTransition,
			ImmediateValueFunction<CashStateMulti, Actions, Demands, Double> immediateValue) {
		this.sampleNum = sampleNum;
		this.distributionsMulti = distributions;
		this.discountFactor = discountFactor;
		this.recursion = recursion;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
	}
	
	public void setSampleNum(int n) {
		this.sampleNum = n;
	}
	
	/**
	 * 
	 * @param iniState
	 * @return simulate sdp in a given number of samples
	 */
	public double simulateSDPGivenSamplNum(CashStateMulti iniState) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; 
			CashStateMulti state = iniState;
			for (int t = 0; t < distributions.length; t++) {
				recursion.getExpectedValue(state);
				Actions actions = new Actions(recursion.getAction(state).getFirstAction(), recursion.getAction(state).getSecondAction());
//				if (t==0)
//					actions = new Actions(14, 16);
				Demands randomDemands = new Demands((int) Math.round(samples[i][t]), (int) Math.round(samples[i][t + distributions.length]));
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, actions, randomDemands);
				state = stateTransition.apply(state, actions, randomDemands);				
			}
			simuValues[i] = sum;			
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}
	
	/**
	 * 
	 * @param iniState
	 * @return simulate sdp in a given number of samples for two product poisson
	 */
	public double simulateSDPGivenSamplNumMulti(CashStateMulti iniState) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamplesMulti(distributionsMulti, sampleNum);
		
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; 
			CashStateMulti state = iniState;
			for (int t = 0; t < distributionsMulti.length; t++) {
				recursion.getExpectedValue(state);
				Actions actions = new Actions(recursion.getAction(state).getFirstAction(), recursion.getAction(state).getSecondAction());
				Demands randomDemands = new Demands((int) Math.round(samples[i][t]), (int) Math.round(samples[i][t + distributionsMulti.length]));
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, actions, randomDemands);
				state = stateTransition.apply(state, actions, randomDemands);				
			}
			simuValues[i] = sum;			
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		System.out.println("\nfinal simulated expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}

}
