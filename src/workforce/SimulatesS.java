package workforce;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.concurrent.CompletionStage;

import sdp.inventory.Recursion;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.BinomialDist;

public class SimulatesS {
	StaffRecursion recursion;
	int T;
	double[] dimissionRate;
	protected StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition; 
	protected ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue; 
	
	public SimulatesS(StaffRecursion recursion, int T, double[] dimissionRate) {
		this.recursion = recursion;
		this.T = T;
		this.stateTransition = recursion.getStateTransitionFunction();
		this.immediateValue = recursion.getImmediateValueFunction();		
		this.dimissionRate = dimissionRate;
	}
	
	public double simulatesS(StaffState iniState, double[][] optimalsS) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		
		int[] sampleNums = new int[T];
		for (int t = 0; t < T; t++) {
			if (t == 0 || t == 1 )
				sampleNums[t] = 10;
			else if (t > 3)
				sampleNums[t] = 1;
			else
				sampleNums[t] = 10;
		}
		
		int sampleTotalNum = Arrays.stream(sampleNums).reduce(1,(a,b) -> a*b);
		int[][] nextInventory = new int[T][];
		double[][] endValues = new double[T][];
		
		for (int t = 0; t < T; t++) {
			int K = sampleNums[t];			
			int sampleNumTot = Arrays.stream(sampleNums).limit(t+1).reduce(1,(a,b) -> a*b);
			endValues[t] = new double[sampleNumTot];
			nextInventory[t] = new int[sampleNumTot];
			int lastStatesLength = t == 0 ? 1 : nextInventory[t-1].length; 
			for (int i = 0; i < lastStatesLength; i++) {
				StaffState thisIniState = t == 0 ? iniState : new StaffState(t + 1, nextInventory[t-1][i]);
				int optQ = thisIniState.iniStaffNum < (int)optimalsS[t][0] ? (int)optimalsS[t][1] - thisIniState.iniStaffNum : 0;
				
				//				int optQ =  Math.max((int)optimalsS[t][1] - thisIniState.iniStaffNum, 0);
				int hireTo = thisIniState.iniStaffNum + optQ;
				int[] randomDemands;
				if (hireTo > 0) {
					BinomialDist dist = new BinomialDist(hireTo, dimissionRate[t]);
					randomDemands = sampling.generateLHSamples(dist, K);
				}
				else {
					randomDemands = new int[K];
					Arrays.fill(randomDemands, 0);
				}
				
				
				for (int j = 0; j < randomDemands.length; j++) {
					try {
						StaffState newState = stateTransition.apply(thisIniState, optQ, randomDemands[j]);
						nextInventory[t][i*K+j] =  newState.iniStaffNum; 
					}
					catch (Exception e) {
						System.out.println(K);
					}
					
					
					if (t > 0)
						endValues[t][i*K+j] = endValues[t-1][i] + immediateValue.apply(thisIniState, optQ, randomDemands[j]);
					else {
						endValues[t][i*K+j] = immediateValue.apply(thisIniState, optQ, randomDemands[j]);
					}
				}				
			}
		}
		
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(endValues[T-1]).sum()/sampleTotalNum;
		System.out.println("\nfinal simulated expected value in " + df2.format(sampleTotalNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}	
}
