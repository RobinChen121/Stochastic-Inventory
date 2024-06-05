package cash.overdraft;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.cash.CashLeadtimeState;
import sdp.cash.CashLeadtimeRecursion;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
*@author: zhenchen
*@date: Mar 1, 2024, 6:21:34 PM
*@desp: 4 periods, demand is 20 each period, running time 50s;
* poisson distribution is faster computed than normal distribution.
* 4 periods, 20 mean demand may be the maximum computational capacity for java, otherwise memory errors;
*
*/

public class SingleProductLeadtime {
	
	public static void main(String[] args) {
		double[] meanDemand = {20, 10, 20, 10};
		
		int T = meanDemand.length;		
		double[] overheadCost = new double[T];
		Arrays.fill(overheadCost, 50);
		double fixOrderCost = 0;
		double variCost = 1;
		double holdingCost = 0;
		double price = 10;
		double salvageValue = 0.5*variCost;
		int leadtime = 1;

		double iniCash = 0;
		double r0 = 0;
		double r1 = 0;
		double r2 = 0.1;
		double r3 = 2; // penalty interest rate for overdraft exceeding the limit
		double limit = 500; // overdraft limit
		double interestFreeAmount = 0;
		double maxOrderQuantity = 30; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 60;
		double minCashState = -200;
		double maxCashState = 300;
		double discountFactor = 1;
		

		// get demand possibilities for each period		
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				//.mapToObj(i -> new NormalDist(meanDemand[i], 0.5*meanDemand[i])) // can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i]))
				//.mapToObj(i -> new GammaDist(meanDemand[0]* beta[1], beta[1]))
				.toArray(Distribution[]::new);
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		// feasible actions
		Function<CashLeadtimeState, double[]> getFeasibleAction = s -> {
			double maxQ = maxOrderQuantity; //- (s.getPeriod()-1)*10; // Math.max(s.iniCash/variCost, 0);//
			if (s.getPeriod() == T)
				maxQ = 0;
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		// when to pay interest is important
		// hold cost is zero
		ImmediateValueFunction<CashLeadtimeState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + state.getPreQ(), randomDemand);
			double variableCost = variCost * action;
			double inventoryLevel = state.getIniInventory() + state.getPreQ() - randomDemand;
			int t = state.getPeriod() - 1;
			double cashBalanceBefore = state.getIniCash()- variableCost - overheadCost[t];// whether plus revenue in this time point
			double interest = 0;
			if (cashBalanceBefore >= 0)
				interest = -r0 * cashBalanceBefore;
			else if(cashBalanceBefore >= -interestFreeAmount)
				interest = 0;
			else if (cashBalanceBefore >= -limit)
				interest = r2 * (-cashBalanceBefore - interestFreeAmount);
			else 
				interest = r3 * (-cashBalanceBefore - limit) + r2 * (limit - interestFreeAmount);

			
			double cashBalanceAfter = cashBalanceBefore - interest + revenue;
			double cashIncrement = cashBalanceAfter - state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			return cashIncrement;
		};
		
		// state transition function
		StateTransitionFunction<CashLeadtimeState, Double, Double, CashLeadtimeState> stateTransition = (state, action,
				randomDemand) -> {
			double nextInventory = Math.max(0, state.getIniInventory() + state.getPreQ() - randomDemand);
			double nextCash = state.getIniCash() + immediateValue.apply(state, action, randomDemand);
			double nextPreQ = action;
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			// cash is integer or not
			//nextCash = Math.round(nextCash * 1) / 1; // affect computation time much
			return new CashLeadtimeState(state.getPeriod() + 1, nextInventory, nextCash, nextPreQ);
		};
		
		/*******************************************************************
		 * Solve
		 */
		CashLeadtimeRecursion recursion = new CashLeadtimeRecursion(pmf, getFeasibleAction, stateTransition, immediateValue);
		int period = 1;
		double iniInventory = 0;
		CashLeadtimeState initialState = new CashLeadtimeState(period, iniInventory, iniCash, 0);
		long currTime = System.currentTimeMillis();
		double opt = recursion.getExpectedValue(initialState);
		System.out.println("final optimal expected value is: " + opt);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		
		double[][] optTable = recursion.getOptTable();	
//		Map<LeadtimeState, Double> cacheValues = recursion.getCacheValues();
		System.out.println("");
	

	}
	
}


