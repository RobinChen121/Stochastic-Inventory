package cash.overdraft;

import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.cash.CashState;
import sdp.cash.CashStateXR;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author chen
 * @email: 15011074486@163.com
 * @Date: 2021 Jan 20 16:53:57
 * @Description: TODO 
 * 
 */
public class CashOverdraftXR {

	public static void main(String[] args) {
		double[] meanDemand = {5, 5, 3};

		double fixOrderCost = 10;
		double variCost = 1;
		double holdingCost = 0;
		double price = 10;
		double salvageValue = 0.5;

		double iniCash = 0;
		double interestRate = 5;
		double minCashRequired = -1000; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 100; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 100;
		double minCashState = -500;
		double maxCashState = 1000;
		double discountFactor = 1;
		

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other distributions
				.toArray(PoissonDist[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		// feasible actions
		Function<CashStateXR, double[]> getFeasibleAction = s -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (s.getIniR() -minCashRequired - fixOrderCost) / variCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};
		
		// immediate value
		ImmediateValueFunction<CashStateXR, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double initCash = state.getIniR() - variCost * state.getIniInventory();
			double cashBalanceBefore = initCash + revenue - fixedCost - variableCost - holdCosts;
			double interest = interestRate * Math.max(-cashBalanceBefore, 0);
			double cashBalanceAfter = cashBalanceBefore - interest;
			double cashIncrement = cashBalanceAfter - state.getIniR();
			return cashIncrement;
		};

		// state transition function
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition = (state, action,
				randomDemand) -> {
			double nextInventory = Math.max(0, state.getIniInventory() + action - randomDemand);
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double holdCosts = holdingCost * Math.max(nextInventory, 0);
			double nextCash = state.getIniCash() + revenue - fixedCost - variableCost - holdCosts;
			nextCash = nextCash - Math.max(-nextCash, 0) * interestRate;
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			// cash is integer or not
			// nextCash = Math.round(nextCash * 100) / 100.00;
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};		

	}

}
