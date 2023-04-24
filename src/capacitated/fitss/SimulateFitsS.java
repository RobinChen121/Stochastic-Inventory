package capacitated.fitss;

import java.util.Arrays;

import sdp.inventory.Recursion;
import sdp.inventory.Simulation;
import sdp.inventory.State;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 13, 2018---2:44:52 PM
*@description:  simulate fitted (s, S) policy for CLSP
*/

public class SimulateFitsS extends Simulation{

	public SimulateFitsS(Distribution[] distributions, int sampleNum, Recursion recursion) {
		super(distributions, sampleNum, recursion);
		// TODO Auto-generated constructor stub
	}

	/**
	 * 
	 * @param iniState
	 * @param optsS
	 * @param maxOrderQuantity
	 * @return simulation results for one level fitted s S policy
	 */
	public double simulateSinglesS(State iniState, double[][] optsS, int maxOrderQuantity){
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; State state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t== 0) 
					optQ = optsS[t][1] - iniState.getIniInventory();
				else 
					optQ =  state.getIniInventory() >= optsS[t][0] ? 0 : Math.min(maxOrderQuantity, optsS[t][1] - state.getIniInventory());;
				sum += immediateValue.apply(state, optQ, samples[i][t]);
				state = stateTransition.apply(state, optQ, samples[i][t]);
			}
			simuValues[i] = sum;
		}
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length;
		System.out.println("\nfinal simulated expected value for fitted one level (s, S) policy is " + simFinalValue);
		return simFinalValue;
	}
	
	/**
	 * 
	 * @param iniState
	 * @param optsS
	 * @param maxOrderQuantity
	 * @return simulation results for two level fitted s S policy
	 */
	public double simulateTwosS(State iniState, double[][] optsS, int maxOrderQuantity){
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] costs = new double[samples.length];
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; State state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t== 0) 
					optQ = optsS[t][1] - iniState.getIniInventory();
				else {
					if (state.getIniInventory() < optsS[t][0])
						optQ = Math.min(maxOrderQuantity, optsS[t][1] - state.getIniInventory());
					else if (optsS[t][0] <= state.getIniInventory() && state.getIniInventory() < optsS[t][2])
						optQ = Math.min(maxOrderQuantity, optsS[t][3] - state.getIniInventory());
					else
						optQ = 0.0;
				}
				sum += immediateValue.apply(state, optQ, samples[i][t]);
				state = stateTransition.apply(state, optQ, samples[i][t]);
				costs[i] = sum;
			}
		}
		double simFinalValue = Arrays.stream(costs).sum()/samples.length;
		System.out.println("\nfinal simulated expected value for fitted two level (s, S) policy is " + simFinalValue);
		return simFinalValue;
	}
	
	/**
	 * 
	 * @param iniState
	 * @param optsS
	 * @param maxOrderQuantity
	 * @return simulation results for three level fitted s S policy
	 */
	public double simulateThreesS(State iniState, double[][] optsS, int maxOrderQuantity){
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] costs = new double[samples.length];
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; State state = iniState;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ;
				if ( t== 0) 
					optQ = optsS[t][1] - iniState.getIniInventory();
				else {
					if (state.getIniInventory() < optsS[t][0])
						optQ = Math.min(maxOrderQuantity, optsS[t][1] - state.getIniInventory());
					else if (optsS[t][0] <= state.getIniInventory() && state.getIniInventory() < optsS[t][2])
						optQ = Math.min(maxOrderQuantity, optsS[t][3] - state.getIniInventory());
					else if (optsS[t][2] <= state.getIniInventory() && state.getIniInventory() < optsS[t][4])
						optQ = Math.min(maxOrderQuantity, optsS[t][5] - state.getIniInventory());
					else
						optQ = 0.0;
				}
				sum += immediateValue.apply(state, optQ, samples[i][t]);
				state = stateTransition.apply(state, optQ, samples[i][t]);
				costs[i] = sum;
			}
		}
		double simFinalValue = Arrays.stream(costs).sum()/samples.length;
		System.out.println("\nfinal simulated expected value for fitted three level (s, S) policy is " + simFinalValue);
		return simFinalValue;
	}
}
