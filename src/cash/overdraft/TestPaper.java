/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: May 16, 2019, 4:30:11 PM
 * @Desc: Test the numerical experiments in paper "Cash-Flow Based Dynamic Inventory Management"
 *
 *
 * 
 */
package cash.overdraft;



import java.util.function.Function;
import java.util.stream.DoubleStream;

import sdp.cash.CashState;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.UniformIntDist;

public class TestPaper {
	
	

	public static void main(String[] args) {
	
		// parameters
		double price = 2000; // selling price 
		double orderCost = 1000; // ordering cost
		double holdCost = 500; // holding cost
		double depositeRate = 0.01; // deposite interest rate
		double loanRate = 0.15; // loan interest rate
		double salvageValue = 600; // salvage value
		
		// demand possible values setting
		double truncationQuantile = 1; 
		int stepSize = 1;
		
		int N = 6; // length of time horizon
		
		double maxOrderQuantity = 200; // maximum ordering quantity in stochastic dynamic programming
		double minInventoryState = 0;
		double maxInventoryState = 300;
		double minCashState = -100; 
		double maxCashState = 1000000;
		
		// distribution: uniform distribution in [0 20]
		Distribution[] distributions = new UniformIntDist[N];
		for (int n = 0; n < N; n++) {
			distributions[n] = new UniformIntDist(0, 20);
		}
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		 
		// feasible actions
		Function<CashState, double[]> getFeasibleAction = p -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (p.getIniCash()) /orderCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};
		
		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double variableCost = orderCost * action;
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdCost * Math.max(inventoryLevel, 0);
			double deposites = depositeRate *  Math.max(state.getIniCash() - variableCost, 0);
			double loanPayed = loanRate * Math.max(variableCost - state.getIniCash(), 0);
			double cashIncrement = revenue - variableCost - holdCosts + deposites - loanPayed;
			double salValue = state.getPeriod() == N ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			return cashIncrement;
		};
	}

}
