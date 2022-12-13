/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 18, 2019, 9:49:16 PM
 * @Desc: This class is to build a stochastic dynamic programming model for a cash constrained problem 
 *        with two products. The main objective of this file is to compute G
 *        
 *
 *
 * 
 */
package cash.multiItem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.DoubleStream;

import sdp.cash.multiItem.RecursionG;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.Recursion;
import sdp.inventory.State;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.inventory.Recursion.OptDirection;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdistmulti.BiNormalDist;


public class MultiItemCashG {


	public static void main(String[] args) {
		double[] price = {2, 10};
		double[] variCost = {1, 2};  // higher margin vs lower margin
		
		double iniCash = 10;  // initial cash
		int[] iniInventory = {0, 0};  // initial inventory
		
		int d = 1; // compute G(d)
		int T = 4;
		
		// mean demand is shape * scale and variance is shape * scale^2
		double[] meanDemands = new double[] {10, 3};
		
		double[][] demand = new double[2][T]; // higher average demand vs lower average demand
		double[] beta = {10, 1}; // higher variance vs lower variance	
		double d1 = meanDemands[0];
		double d2 = meanDemands[1];
		for (int t = 0; t < T; t++) {
			demand[0][t] = d1;
			demand[1][t] = d2;
		}
		
		
		double[] salPrice = Arrays.stream(variCost).map(a -> a*0.5).toArray();
		
		double truncationQuantile = 0.9999;
		int stepSize = 1;			
		int maxInventoryState = 200;
		int Qbound = 40;
		double discountFactor = 1;
		

		
		// get shape possibilities for a product in each period
		GammaDist[] distributions =  new GammaDist[T]; // normal dist for one product
		for (int t = 0; t < T; t++)
			distributions[t] = new GammaDist(demand[d-1][t]* beta[d-1], beta[d-1]);

		
		// build action list for this item
		Function<State, double[]> buildActionList = s -> {
			return DoubleStream.iterate(0, i -> i + stepSize).limit(Qbound + 1).toArray();
		};
		
		// Immediate Value Function	      
		ImmediateValueFunction<State, Double, Double, Double> immediateValue
		= (IniState, action, randomDemand) -> {
			double revenue = 0;
			revenue = (price[d - 1] - variCost[d -1]) * Math.min(IniState.getIniInventory() + action, randomDemand);
			if (IniState.getPeriod() == T) {
				revenue += (salPrice[d - 1] - variCost[d -1]) * Math.max(IniState.getIniInventory() + action - randomDemand, 0);
			}
			return revenue;
		};
	    	
		// State Transition Function
		// need change
		StateTransitionFunction<State, Double, Double, State> stateTransition = (IniState, action, randomDemand) -> {
			double endInventory = IniState.getIniInventory() + action - randomDemand;
			endInventory = Math.max(0, endInventory);
			endInventory = endInventory > maxInventoryState ? maxInventoryState : endInventory;           
			endInventory = (int) endInventory;
			return new State(IniState.getPeriod() + 1, endInventory);
		};
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		
		/*******************************************************************
		 * Solve
		 */
		RecursionG recursion = new RecursionG(pmf, buildActionList,
				                             stateTransition, immediateValue);
		int period = 1;
		State iniState = new State(period, iniInventory[d - 1]);
		long currTime = System.currentTimeMillis();
		double finalValue = iniCash + recursion.getExpectedValue(iniState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first priod is :  Q = " + recursion.getAction(iniState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		System.out.println("a* in each period:");
		double[] optY = recursion.getOptY();
		optY[0] = iniState.getIniInventory() + recursion.getAction(iniState);
		System.out.println(Arrays.toString(optY));
		
		
		

		
		
	}

}
