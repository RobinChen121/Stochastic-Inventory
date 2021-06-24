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

public class CashSimulationMultiXR {
	int sampleNum;

	
	Distribution[][] distributionsMulti;
	
	double discountFactor;
	
	CashRecursionMultiXR recursion;
	
	StateTransitionFunction<CashStateMultiXR, double[], double[], CashStateMultiXR> stateTransition;
	ImmediateValueFunction<CashStateMultiXR, double[], double[], Double> immediateValue;
	

	
	public CashSimulationMultiXR(int sampleNum, Distribution[][] distributions, double discountFactor, 
			CashRecursionMultiXR recursion, StateTransitionFunction<CashStateMultiXR, double[], double[], CashStateMultiXR> stateTransition,
			ImmediateValueFunction<CashStateMultiXR, double[], double[], Double> immediateValue) {
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
	public double simulateSDPGivenSamplNum(CashStateMultiXR iniState) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributionsMulti, sampleNum);
		
//		double sumD = 0;
//		for(int i = 0; i < samples.length; i++) {
//			sumD += samples[i][3];
//		}
//		double meanD = sumD/samples.length;
		
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; 
			CashStateMultiXR state = iniState;
			for (int t = 0; t < distributionsMulti[0].length; t++) {
				recursion.getExpectedValue(state);
				double[] actions = new double[] {recursion.getAction(state)[0], recursion.getAction(state)[1]};
				double[] randomDemands = new double[] {samples[i][t* 2], samples[i][t* 2 + 1]};
				sum += Math.pow(discountFactor, t) * immediateValue.apply(state, actions, randomDemands);
				state = stateTransition.apply(state, actions, randomDemands);				
			}
			simuValues[i] = sum;			
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniR;
		System.out.println("\nfinal simulated expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}
	


}
