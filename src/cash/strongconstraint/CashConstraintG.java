
package cash.strongconstraint;

import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.IntStream;

import sdp.cash.RecursionG;
import sdp.cash.StateY;
import sdp.inventory.GetPmf;
import sdp.inventory.State;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Apr 17, 2020, 4:37:46 PM
 * @Desc: 
 *Cash constraint function G() for Chao's paper
 *
 * 
 */

public class CashConstraintG {
	public static void main(String[] args) {
		double[] meanDemand = {8, 8, 8};

		double variCost = 5;
		double price = 10;
		double depositeRate = 0;
		double salvageValue = 1;
		double truncationQuantile = 0.99;
		int stepSize = 1;

		double coe = 1;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i]))			
				//.mapToObj(i -> new NormalDist(meanDemand[i], coe * Math.sqrt(meanDemand[i]))) // can be changed to other distributions				
				//.mapToObj(i -> new GammaDist(meanDemand[i], 2))
				.toArray(Distribution[]::new);
				
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		/*******************************************************************
		 * Solve
		 */

		RecursionG recursion = new RecursionG(pmf, distributions, 
				price, variCost, depositeRate, salvageValue);
		long currTime = System.currentTimeMillis();	
		double[] optY = recursion.getOptY();
		System.out.println("a* in each period: ");
		System.out.println(Arrays.toString(optY));
		double time = (System.currentTimeMillis() - currTime) / 1000.0; // if 1000, then it will be integer
		System.out.printf("running time is %.3f s", time);

	}

}
