package cash.overdraft;

import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.cash.CashLeadtimeState;
import sdp.cash.CashLeadtimeRecursion;
import sdp.inventory.GetPmf;
import sdp.inventory.LeadtimeRecursion;
import sdp.inventory.LeadtimeState;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
*@author: zhenchen
*@date: Mar 1, 2024, 6:21:34 PM
*@desp: TODO
*
*/

public class SingleProductLeadtime {
	
	public static void main(String[] args) {
		double[] meanDemand = {20, 35, 20};
		
		double[] overheadCost = {100, 100, 100};
		double fixOrderCost = 0;
		double variCost = 1;
		double holdingCost = 0;
		double price = 10;
		double salvageValue = 0;
		int leadtime = 1;

		double iniCash = 0;
		double r0 = 0.01;
		double r1 = 0;
		double r2 = 0.1;
		double r3 = 1; // penalty interest rate for overdraft exceeding the limit
		double limit = 1000; // overdraft limit
		double interestFreeAmount = 0;
		double maxOrderQuantity = 100; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 100;
		double minCashState = -200;
		double maxCashState = 800;
		double discountFactor = 1;
		

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other distributions
				.toArray(PoissonDist[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		// feasible actions
		Function<CashLeadtimeState, double[]> getFeasibleAction = s -> {
			double maxQ = maxOrderQuantity; // Math.max(s.iniCash/variCost, 0);//
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
			//nextCash = Math.round(nextCash * 10) / 10; // affect computation time much
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


