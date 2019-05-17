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



import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.cash.CashRecursion;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.probdist.UniformIntDist;

public class TestPaper {
	
	

	public static void main(String[] args) {
			
		double iniInventory = 0;
		int iniCash = 0;
		
		// parameters
		double price = 2000; // selling price 
		double orderCost = 1000; // ordering cost
		double holdCost = 500; // holding cost
		double depositeRate = 0.01; // deposite interest rate
		double loanRate = 0.15; // loan interest rate
		double salvageValue = 400; // salvage value
		
		// demand possible values setting
		double truncationQuantile = 1; 
		int stepSize = 1;
		
		int N = 12; // length of time horizon
		
		double maxOrderQuantity = 100; // maximum ordering quantity in stochastic dynamic programming
		double minInventoryState = 0;
		double maxInventoryState = 300;
		double minCashState = -1000000; 
		double maxCashState = 1000000;
		double minCashRequired = -1000000;
		
		double discountFactor = 1;
		
		// distribution: uniform distribution in [0 20]
		Distribution[] distributions = new UniformIntDist[N];
		for (int n = 0; n < N; n++) {
			distributions[n] = new UniformIntDist(0, 20);
		}
		
//		double[] meanDemand = {10, 7, 20, 5, 15};
//		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(N)
//				//.mapToObj(i -> new NormalDist(meanDemand[i], 0.25 * meanDemand[i])) // can be changed to other distributions
//				.mapToObj(i -> new PoissonDist(meanDemand[i]))
//				.toArray(Distribution[]::new);
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		 
		// feasible actions
		Function<CashState, double[]> getFeasibleAction = p -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (p.getIniCash() - minCashRequired) /orderCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};
		
		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double variableCost = orderCost* action;
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = state.getPeriod() == N ? 0 : holdCost * Math.max(inventoryLevel, 0);
			double deposites = depositeRate *  Math.max(state.getIniCash() - variableCost, 0);
			double loanPayed = loanRate * Math.max(variableCost - state.getIniCash(), 0);
			double cashIncrement = revenue - variableCost - holdCosts + deposites - loanPayed;
			double salValue = state.getPeriod() == N ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			return cashIncrement;
		};
		
		// state transition function
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition = (state, action,
				randomDemand) -> {
					double nextInventory = Math.max(0, state.getIniInventory() + action - randomDemand);
					double nextCash = state.getIniCash() + immediateValue.apply(state, action, randomDemand);
					nextCash = nextCash > maxCashState ? maxCashState : nextCash;
					nextCash = nextCash < minCashState ? minCashState : nextCash;
					nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
					nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
					// cash is integer or not
					
					if (state.getPeriod() > 2)
						nextCash = Math.round(nextCash * 0.0001) / 0.0001; 
					return new CashState(state.getPeriod() + 1, nextInventory, (int) nextCash);
		};
		
		/*******************************************************************
		 * Solve
		 */
		CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition,
				immediateValue, discountFactor);
		int period = 1;		
		CashState initialState = new CashState(period, iniInventory, iniCash);
		long currTime = System.currentTimeMillis();
		recursion.setTreeMapCacheAction();
		double finalValue = iniCash + recursion.getExpectedValue(initialState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		double[][] optTable = recursion.getOptTable();
		
		// compute alpha, beta
		double[] a = new double[N];
		double[] b = new double[N];
		double[] alpha = new double[N];
		double[] beta = new double[N];
		for (int i = 0; i < N; i++) {
			if (i == N - 1) {
				a[i] = (price - orderCost * (1 + loanRate)) / (price - salvageValue);
				b[i] = (price - orderCost * (1 + depositeRate)) / (price - salvageValue);
			}
			else {
				a[i] = (price - orderCost * (1 + loanRate)) / (price + holdCost);
				b[i] = (price - orderCost * (1 + depositeRate)) / (price + holdCost);
			}
			alpha[i] = distributions[i].inverseF(a[i]);
			beta[i] = distributions[i].inverseF(b[i]);
		}
		System.out.println("alpha: " + Arrays.toString(alpha));
		System.out.println("beta: " + Arrays.toString(beta));
		
	}

}
