package milp;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

import sdp.cash.CashState;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import umontreal.ssj.functions.ShiftedMathFunction;
import umontreal.ssj.probdist.Distribution;

public class SimulateChanceCash {
	double iniCash;
	double iniI;
	double[] price;
	double variCostUnit;
	double salvageValueUnit;
	double[] overheadCost;
	double serviceRate;
	protected StateTransitionFunction<CashState, Double, Double, CashState> stateTransition; 
	protected ImmediateValueFunction<CashState, Double, Double, Double> immediateValue; 
	double discountFactor;
	Distribution[] distributions;
	int[] sampleNums;
	double[][] scenarios;
	double holdCostUnit;
	
	public SimulateChanceCash(Distribution[] distributions, double iniCash, double iniI, double[] price, double variCostUnit, double salvageValueUnit,
			double holdCostUnit, double[] overheadCost, double serviceRate, StateTransitionFunction<CashState, Double, Double, CashState> stateTransition,
			ImmediateValueFunction<CashState, Double, Double, Double> immediateValue, double discountFactor,
			int[] sampleNums, double[][] scenarios) {
		this.iniCash = iniCash;
		this.iniI = iniI;
		this.price = price;
		this.variCostUnit = variCostUnit;
		this.salvageValueUnit = salvageValueUnit;
		this.overheadCost = overheadCost;
		this.serviceRate = serviceRate;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		this.discountFactor = discountFactor;		
		this.distributions = distributions;
		this.sampleNums = sampleNums;
		this.scenarios = scenarios;
		this.holdCostUnit = holdCostUnit;
	}
	
	public double simulateSDPGivenSamplNum(CashState iniState, double iniQ, int sampleNum) {
		Sampling.resetStartStream();
		Sampling sampling = new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];	
		
		for (int i = 0; i < samples.length; i++) {
			double sum = 0;
			CashState state = iniState;
			int T = sampleNums.length;
			for (int t = 0; t < samples[0].length; t++)
			{
				if (t == 0) {
					double optQ = iniQ;
					double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
					double thisValue = immediateValue.apply(state, optQ, randomDemand);
					sum += Math.pow(discountFactor, t) * thisValue;
					state = stateTransition.apply(state, optQ, randomDemand);
				}
				else {
					double tIniCash = state.iniCash;
					double tIniI = state.getIniInventory();
					PositiveCashChance model = new PositiveCashChance(Arrays.copyOfRange(distributions, t, T), Arrays.copyOfRange(sampleNums, t, T), tIniCash, tIniI, price, variCostUnit, 
							salvageValueUnit, holdCostUnit, overheadCost, serviceRate, Arrays.copyOfRange(scenarios, t, T));
					
					double[] result = new double[3];
					try {
						result = model.solveSort(); // sort or not sort
					} catch (Exception e) {
						System.out.println(result);
					}
					
					double optQ = result[0];
					double randomDemand = Math.round(samples[i][t]); // integer samples to test sdp
					double thisValue = immediateValue.apply(state, optQ, randomDemand);
					sum += Math.pow(discountFactor, t) * thisValue;
					state = stateTransition.apply(state, optQ, randomDemand);
				}
				simuValues[i] = sum;
			}
		}
		//DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length + iniState.iniCash;
		//System.out.println("\nfinal simulated expected value for chanced SAA in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}
	
	


}
