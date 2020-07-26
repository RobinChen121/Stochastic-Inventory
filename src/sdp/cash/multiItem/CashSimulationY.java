/**
 * @date: Jul 19, 2020
 */
package sdp.cash.multiItem;

import java.text.DecimalFormat;
import java.util.Arrays;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.inventory.StateTransition.StateTransitionFunctionV;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 19, 2020
 * @Desc: simulate the ordering policy of cash constrained two product problem by Younes
 *
 */
public class CashSimulationY {
	int sampleNum;	
	Distribution[][] distributionsMulti;	
	double discountFactor;	
	CashRecursionV recursion;
	
	StateTransitionFunctionV<CashStateMultiYR, double[], CashStateMulti> stateTransition;
	
	public CashSimulationY(int sampleNum, Distribution[][] distributions, double discountFactor, 
			CashRecursionV recursion, StateTransitionFunctionV<CashStateMultiYR, double[], CashStateMulti> stateTransition) {
		this.sampleNum = sampleNum;
		this.distributionsMulti = distributions;
		this.discountFactor = discountFactor;
		this.recursion = recursion;
		this.stateTransition = stateTransition;
	}
	
	

	public void setSampleNum(int n) {
		this.sampleNum = n;
	}
	
	public double simulateSDPGivenSamplNum(CashStateMulti iniState, double[] variCost) {
		Sampling.resetStartStream();
		double[][] samples = Sampling.generateLHSamples(distributionsMulti, sampleNum);
		
		double[] simuValues = new double[samples.length];		
		for (int i = 0; i < samples.length; i++) {
			double finalValue = 0; 
			CashStateMulti state = iniState;
			for (int t = 0; t < distributionsMulti[0].length; t++) {
				CashStateR stateR = new CashStateR(t + 1, state.getIniCash() + variCost[0] * state.getIniInventory1() + variCost[1] * state.getIniInventory2());
				double[] actionY = recursion.getYStar(stateR);
				double[] actions = new double[] {0, 0};
				double[] randomDemands = new double[] {samples[i][t* 2], samples[i][t* 2 + 1]};
				double alpha = 0;
				if (state.getIniInventory1() < actionY[0]+0.1 && state.getIniInventory2() < actionY[1]
						&& variCost[0] * actionY[0] + variCost[1] * actionY[1] < stateR.iniR+0.1) {
					actions[0] = actionY[0]; actions[1] = actionY[1];
				}
				else if (state.getIniInventory1() > actionY[0]-0.1 && state.getIniInventory2() > actionY[1]-0.1) {
					actions[0] = state.getIniInventory1(); actions[1] = state.getIniInventory2();
				}
				else if (state.getIniInventory1() > actionY[0]-0.1 && state.getIniInventory2() < actionY[1]+0.1) {
					double x1 = state.getIniInventory1();
					actions[0] = state.getIniInventory1(); actions[1] = Math.min(actionY[1], (stateR.iniR - x1*variCost[0]) / variCost[1]);
				}
				else if (state.getIniInventory1() < actionY[0]+0.1 && state.getIniInventory2() > actionY[1]-0.1){
					double x2 = state.getIniInventory2();
					actions[0] = Math.min(actionY[0], (stateR.iniR -  x2*variCost[1]) / variCost[0]); actions[1] = state.getIniInventory2();
				}
				else if (state.getIniInventory1() < actionY[0]+0.1 && state.getIniInventory2() < actionY[1]+0.1
						&& variCost[0] * actionY[0] + variCost[1] * actionY[1] > stateR.iniR-0.1) {
					alpha =  recursion.getAlpha(stateR);
					actions[0] = alpha * stateR.iniR / variCost[0]; 
					actions[1] = (1 - alpha) * stateR.iniR / variCost[1];				
				}
				
				CashStateMultiYR thisState = new CashStateMultiYR(t + 1, actions[0], actions[1], stateR.iniR);
				CashStateMulti newState = stateTransition.apply(thisState, randomDemands);					
				state = newState;
				if (t == distributionsMulti[0].length - 1)
					finalValue = recursion.boundFinalCash.apply(newState);
			}
			simuValues[i] = finalValue;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length;
		System.out.println("\nfinal simulated expected value for this policy in " + df2.format(sampleNum) + " samples is: " + df2.format(simFinalValue));
		return simFinalValue;
	}

}
