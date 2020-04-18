/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Apr 17, 2020, 4:37:46 PM
 * @Desc: 
 *Cash constraint function G() for Chao's paper
 *
 * 
 */
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
import umontreal.ssj.probdist.NormalDist;

public class CashConstraintG {
	public static void main(String[] args) {
		double[] meanDemand = {10, 10, 10, 10};
		double iniInventory = 0;
		double iniCash = 7;
		double fixOrderCost = 0;
		double variCost = 1;
		double price = 1.3;
		double depositeRate = 0.1;
		double salvageValue = 0.5;
		double holdingCost = 0;	

		double overheadCost = 0; // costs like wages or rents which is required to pay in each period
		double overheadRate = 0; // rate from revenue to pay overhead wages
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.99;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;	
		double discountFactor = 1 + depositeRate;

		double coe = 1;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new NormalDist(meanDemand[i], coe * Math.sqrt(meanDemand[i]))) // can be changed to other distributions
				//.mapToObj(i -> new PoissonDist(meanDemand[i]))
				.toArray(Distribution[]::new);
		
		double aNStar = distributions[T - 1].inverseF((price - (1 + depositeRate) * variCost) / (price - salvageValue)); // a*_N
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		// feasible actions
		Function<StateY, double[]> getFeasibleAction = s -> {
			double[] actionY = new double[(int)maxOrderQuantity + 1];
			for (int i = 0; i < (int) maxOrderQuantity; i++)
				actionY[i] = i;
			return actionY;
		};

		// Immediate Value Function	   
		BiFunction<StateY, Double, Double> immediateValue
		= (IniState, randomDemand) -> {
			double revenue = 0;
			revenue = (price - variCost) * Math.min(IniState.getIniY(), randomDemand);
			double cost = depositeRate * variCost * IniState.getIniY();
			if (IniState.getPeriod() == T)
				revenue += (salvageValue - variCost) * Math.max(aNStar, Math.max(IniState.getIniY() - randomDemand, 0));
			return revenue - cost;
		};		


		StateTransitionFunction<StateY, Double, Double, StateY> stateTransition = (IniState, nextOptY, randomDemand) -> {
			double endY = Math.max(0, IniState.getIniY() - randomDemand);
			return new StateY(IniState.getPeriod() + 1, Math.max(nextOptY, endY));
		};


		/*******************************************************************
		 * Solve
		 */

		RecursionG recursion = new RecursionG(pmf, aNStar, getFeasibleAction,
                stateTransition, immediateValue);
		int period = 1;
		StateY iniState = new StateY(period, 0); // random y = 0
		recursion.getExpectedValue(iniState);
		long currTime = System.currentTimeMillis();
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		double[] optY = recursion.getOptY();
		System.out.println("a* in each period: ");
		System.out.println(Arrays.toString(optY));



	}

}
