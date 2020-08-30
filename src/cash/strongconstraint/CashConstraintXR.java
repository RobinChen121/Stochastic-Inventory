package cash.strongconstraint;

import java.awt.print.Printable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import cash.strongconstraint.FindsCS.FindCCrieria;
import sdp.cash.CashRecursionXR.OptDirection;
import sdp.cash.CashRecursionXR;
import sdp.cash.CashSimulationXR;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToCsv;
import sdp.write.WriteToExcel;
import sdp.cash.CashStateXR;
import sdp.cash.RecursionG;
import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date 2020, Feb 22th, 12:31:10 pm
 * @Description a multi-period stochastic inventory problem to implement the method in Chao's paper:
 *              Dynamic inventory management with cash flow constraints (2008) in Naval Research Logistics.
 *              There is no fixed ordering cost and it is a multi-period newsvendor problem.
 *              
 *              In this dynamic programming, states are inventory level x and working capital R = w+cx, 
 *              action is order-up-to level y.
 *
 */

public class CashConstraintXR {
	
	public static void main(String[] args) {
		double[] meanDemand = {8, 8, 8, 8};
		double iniInventory = 0;
		double iniCash = 30;
		double fixOrderCost = 0;
		double variCost = 2;
		double price = 4;
		double depositeRate = 0;
		double salvageValue = 1;
		double holdingCost = 0;	
		FindCCrieria criteria = FindCCrieria.XRELATE;		
		double overheadCost = 0; // costs like wages or rents which is required to pay in each period
		double overheadRate = 0; // rate from revenue to pay overhead wages
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash

		double truncationQuantile = 0.99;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;
		
		double discountFactor = 1;

		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				//.mapToObj(i -> new NormalDist(meanDemand[i], Math.sqrt(meanDemand[i]))) // can be changed to other distributions
				//.mapToObj(i -> new PoissonDist(meanDemand[i]))
				.mapToObj(i -> new GammaDist(meanDemand[i], 2))
				.toArray(Distribution[]::new);
	
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
			

		// feasible actions
		Function<CashStateXR, double[]> getFeasibleAction = s -> {
			double maxY = s.getIniR() / variCost < s.getIniInventory() ? s.getIniInventory() : s.getIniR() / variCost;
			int length = (int) (maxY - s.getIniInventory()) + 1;
			return DoubleStream.iterate(s.getIniInventory(), i -> i + stepSize).limit(length).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashStateXR, Double, Double, Double> immediateValue = (state, actionY, randomDemand) -> {
			double revenue = price * Math.min(actionY, randomDemand);
			double action = actionY - state.getIniInventory();
			double fixedCost = actionY > state.getIniInventory() ? fixOrderCost : 0;
			double variableCost = variCost * action;
			double initCash = state.getIniR() - variCost * state.getIniInventory();
			double deposite = (initCash- fixedCost - variableCost) * (1 + depositeRate); // (1+d)(S-cy)
			double inventoryLevel = actionY - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);			
			double cashIncrement = (1 - overheadRate)*revenue + deposite - holdCosts - overheadCost - initCash;
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			return cashIncrement;
		};

		// state transition function
		StateTransitionFunction<CashStateXR, Double, Double, CashStateXR> stateTransition = (state, actionY,
				randomDemand) -> {
			if (randomDemand < 0)
				System.out.println(randomDemand);
			double nextInventory = Math.max(0, actionY - randomDemand);
			double initCash = state.getIniR() - variCost * state.getIniInventory();
			double nextCash = initCash + immediateValue.apply(state, actionY, randomDemand);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			// cash is integer or not
			nextCash = Math.round(nextCash * 1) / 1; 
			double nextR = nextCash + variCost * nextInventory;
			
			return new CashStateXR(state.getPeriod() + 1, nextInventory, nextR, variCost);
		};

		/*******************************************************************
		 * Solve
		 */
		CashRecursionXR recursion = new CashRecursionXR(OptDirection.MAX, pmf, getFeasibleAction, stateTransition,
				immediateValue, discountFactor);
		int period = 1;		
		CashStateXR initialState = new CashStateXR(period, iniInventory, iniCash, variCost);
		long currTime = System.currentTimeMillis();
		recursion.setTreeMapCacheAction();
		double finalValue = iniCash + recursion.getExpectedValue(initialState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000.0;
		System.out.println("running time is " + time + "s");

		
		/*******************************************************************
		 * Simulating sdp results
		 * parameter vales like price, variCost, holdingCost etc.
		 * are only for compute L(y), not very necessary
		 */
		int sampleNum = 10000;
		
		CashSimulationXR simuation = new CashSimulationXR(distributions, sampleNum, recursion, discountFactor, 
				fixOrderCost, price, variCost, holdingCost, salvageValue); // no need to add overheadCost in this class
		double simFinalValue = simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001; 
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);
		
		/*******************************************************************
		 * get optimal table of SDP, 
		 * and output it to a excel file
		 */
//		System.out.println("");
//		double[][] optTable = recursion.getOptTable();
//		WriteToExcel wr = new WriteToExcel();
//		String headString =  "period" + "\t" + "x" + "\t" + "S" + "\t" + "R" + "\t" + "y";
//		wr.writeArrayToExcel(optTable, "optTable.xls", headString);
		
		
		/*******************************************************************
		 * get a* in each period
		 */
		RecursionG recursion2 = new RecursionG(pmf, distributions, 
				price, variCost, depositeRate, salvageValue);
		currTime = System.currentTimeMillis();	
		double[] optY = recursion2.getOptY();
		System.out.println();
		time = (System.currentTimeMillis() - currTime) / 1000.0; // if 1000, then it will be integer
		System.out.println("a* in each period: ");
		System.out.println(Arrays.toString(optY));
		System.out.printf("running time is %.3f s", time);
		System.out.println();
		
		/*******************************************************************
		 * simulate a* in each period
		 */
		simuation.simulateAStar(optY, initialState);
		
		
	}			
}
