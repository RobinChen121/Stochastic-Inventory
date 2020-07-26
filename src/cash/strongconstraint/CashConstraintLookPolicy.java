package cash.strongconstraint;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import cash.strongconstraint.FindsCS.FindCCrieria;
import milp.MipCashConstraint;
import sdp.cash.CashRecursion;
import sdp.cash.CashRecursion.OptDirection;
import sdp.cash.CashSimulation;
import sdp.inventory.CheckKConvexity;
import sdp.inventory.Drawing;
import sdp.inventory.GetPmf;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToCsv;
import sdp.write.WriteToExcel;
import sdp.cash.CashState;
import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date 2018, March 3th, 6:31:10 pm
 * @Description find the characteristics of the optimal policy;
 * 
 *
 */

public class CashConstraintLookPolicy {
	

	public static void main(String[] args) {
		double[] meanDemand = {2, 3, 8};
		//double[] meanDemand = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
		double iniInventory = 8;
		double iniCash = 100;
		double fixOrderCost = 10;
		double variCost = 1;
		double price = 8;
		double depositeRate = 0;
		double salvageValue = 0.5;
		double holdingCost = 0;	
		FindCCrieria criteria = FindCCrieria.XRELATE;		
		double overheadCost = 0; // costs like wages or rents which is required to pay in each period
		double overheadRate = 0; // rate from revenue to pay overhead wages
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;
		
		double discountFactor = 1;
		
		double xmin = 0; double xmax = 10;
		double Rmin = fixOrderCost; double Rmax = 40;
		int row = 0;
		int column = 0;
		int columnNum  = (int) (Rmax - Rmin + 1) + 1;
		int rowNum= (int) ((xmax - xmin + 1)/1) + 1; // ((Rmax - Rmin + 1)/2) + 2;
		double[][] resultTable = new double[rowNum][columnNum];
		
		for (iniInventory = xmax; iniInventory >= xmin; iniInventory--) {
			column = 0;
			for (iniCash = Rmin; iniCash <= Rmax; iniCash = iniCash + 1) {
	
		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				//.mapToObj(i -> new NormalDist(meanDemand[i], Math.sqrt(meanDemand[i]))) // can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i]))
				.toArray(Distribution[]::new);
	
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
			

		// feasible actions
		Function<CashState, double[]> getFeasibleAction = s -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (s.getIniCash() - overheadCost - fixOrderCost) / variCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double deposite = (state.getIniCash() - fixedCost - variableCost) * (1 + depositeRate);
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = (1 - overheadRate)*revenue + deposite - holdCosts - overheadCost - state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
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
			nextCash = Math.round(nextCash * 1) / 1; 
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
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
		double optQ = recursion.getAction(initialState);
		System.out.println("optimal order quantity in the first priod is : " + optQ);
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");

		System.out.println("initial inventory is " + iniInventory);
		System.out.println("initial cash is " + iniCash);
		resultTable[0][column + 1] = iniCash;
		resultTable[row + 1][0] = iniInventory;
		resultTable[row + 1][column + 1] = optQ; // finalValue - iniCash;
		System.out.println("**********************************************************");


		
//		/*******************************************************************
//		 * Simulating sdp results
//		 * parameter vales like price, variCost, holdingCost etc.
//		 * are only for compute L(y), not very necessary
//		 */
//		int sampleNum = 10000;
//		
//		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor, 
//				fixOrderCost, price, variCost, holdingCost, salvageValue); // no need to add overheadCost in this class
//		double simFinalValue = simuation.simulateSDPGivenSamplNum(initialState);
//		double error = 0.0001; 
//		double confidence = 0.95;
//		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);
//				
			column++;
			}
		row++;
		}
		WriteToExcel wr = new WriteToExcel();
		wr.writeArrayToExcel(resultTable, "resultTable.xls");
		
	}			
}
