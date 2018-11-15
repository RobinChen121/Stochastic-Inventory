package milp;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.stream.IntStream;

import javax.naming.InitialContext;

import sdp.inventory.State;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;

/**
* @author Zhen Chen
* @date: 2018年11月14日 下午5:38:17  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  this is a class to simulate MILP results
*/

public class Simulation {
	double[] meanDemand; 
	double[] sigma;
	double iniInventory;	
	double fixOrderCost;
	double variCost;
	double holdingCost;
	double penaltyCost;
	
	public Simulation(double[] meanDemand, double[] sigma, double iniInventory, Double fixOrderCost, double variCost, double holdingCost,
			double penaltyCost) {
		this.meanDemand = meanDemand;
		this.sigma = sigma;
		this.iniInventory = iniInventory;
		this.fixOrderCost = fixOrderCost;
		this.variCost = variCost;
		this.holdingCost = holdingCost;
		this.penaltyCost = penaltyCost;		
	}
	
	public double simulatesS(double[][] sS, int sampleNum) {
		Sampling.resetStartStream();
		
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new NormalDist(meanDemand[i], 0.25*meanDemand[i])) // can be changed to other distributions
				.toArray(Distribution[]::new);

		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
		double[] simuValues = new double[samples.length];	
		for (int i = 0; i < samples.length; i++) {
			double sum = 0; 
			double initialInventory = iniInventory;
			for (int t = 0; t < samples[0].length; t++)
			{
				double optQ = initialInventory >= sS[t][0] ? 0 : sS[t][1] - initialInventory;
				double randomDemand = samples[i][t]; // integer samples to test sdp
				double fixCost = optQ > 0 ? fixOrderCost : 0;
				double variCosts = variCost * optQ;
				double inventory = initialInventory + optQ - randomDemand;
				double holdCosts = holdingCost * Math.max(0, inventory);
				double penaCosts = penaltyCost * Math.max(0, -inventory);
				sum = sum + fixCost + variCosts + holdCosts + penaCosts;
				initialInventory = inventory;
			}
			simuValues[i] = sum;
		}
		DecimalFormat df2 = new DecimalFormat("###,###");
		double simFinalValue = Arrays.stream(simuValues).sum()/samples.length;
		System.out.println("final simulated expected value in " + df2.format(sampleNum) + " samples is: " + simFinalValue);
		return simFinalValue;
	}

}


