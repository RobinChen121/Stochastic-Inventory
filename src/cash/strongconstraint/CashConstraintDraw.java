package cash.strongconstraint;

import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import cash.strongconstraint.FindsCS.FindCCrieria;
import sdp.cash.CashRecursion;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.CheckKConvexity;
import sdp.inventory.Drawing;
import sdp.inventory.GetPmf;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToExcel;
import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;


/**
 * @author chen zhen
 * @version 2018, April 15th, 10:32:46 am
 * @Description: drawing pictures, see K convexity;
 * an example of two local minimums:
 * 
double[] meanDemand = { 6.6, 2, 21.8};
		double iniCash = 33;
		double iniInventory = 0;
		double fixOrderCost = 20;
		double variCost = 1;
		double price = 4;
		double salvageValue = 0;
		FindCCrieria criteria = FindCCrieria.XRELATE;
		double holdingCost = 1;	
 *  
 * 
 *  whether order or not can be found by the relation of GA(0) and GB(y*)- K
 *  
 */
public class CashConstraintDraw {


	public static void main(String[] args) {
		double[] meanDemand = {7, 2, 6};
		//double[] meanDemand = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
		double iniCash = 16;
		double iniInventory = 0;
		double fixOrderCost = 15;
		double variCost = 1;
		double price = 6;
		double salvageValue = 0;
		double holdingCost = 0;	
		FindCCrieria criteria = FindCCrieria.XRELATE;		
		double overheadCost = 0; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash
		

		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 500;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 2000;	
		double discountFactor = 1;
		
		boolean isForDrawGy = true;
		// get demand possibilities for each period
		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i])) // can be changed to other distributions
				.toArray(PoissonDist[]::new);
		
//		double[] values = {6, 7};
//		double[] probs = {0.95, 0.05};
//		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
//		.mapToObj(i -> new DiscreteDistribution(values, probs, values.length)) // can be changed to other distributions
//		.toArray(DiscreteDistribution[]::new);	
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		// feasible actions
		Function<CashState, double[]> getFeasibleAction = s -> {
			double maxQ = (int) Math.min(maxOrderQuantity,
					Math.max(0, (s.getIniCash() - overheadCost - fixOrderCost) / variCost));
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};

		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			double revenue = 0;
			double fixedCost = 0;
			double variableCost = 0;
			double inventoryLevel = 0;
			revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
			fixedCost = action > 0 ? fixOrderCost : 0;
			variableCost = variCost * action;
			inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue - fixedCost - variableCost - holdCosts;
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			return cashIncrement;
		};

		// state transition function
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition = (state, action,
				randomDemand) -> {
			double nextInventory = state.getIniInventory() + action - randomDemand;
			double nextCash = state.getIniCash() + immediateValue.apply(state, action, randomDemand);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			// cash is integer or not
			// nextCash = Math.round(nextCash * 100) / 100.00;
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
		double finalCash = iniCash + recursion.getExpectedValue(initialState);
		System.out.println("final optimal cash is: " + finalCash);
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");

		/*******************************************************************
		 * Simulating sdp results
		 */
		int sampleNum = 10000;
		CashSimulation simuation = new CashSimulation(distributions, sampleNum, recursion, discountFactor, fixOrderCost, price, variCost, holdingCost, salvageValue);
		double simFinalValue = simuation.simulateSDPGivenSamplNum(initialState);
		double error = 0.0001;
		double confidence = 0.95;
		simuation.simulateSDPwithErrorConfidence(initialState, error, confidence);
		
		/*******************************************************************
		 * Find (s, C, S) and simulate
		 */
		System.out.println("");
		double[][] optTable = recursion.getOptTable();
		FindsCS findsCS = new FindsCS(iniCash, distributions, fixOrderCost, price, variCost, holdingCost, salvageValue);
		double[][] optsCS = findsCS.getsC12S(optTable, overheadCost, criteria);
		Map<State, Double> cacheC1Values = new TreeMap<>();
		Map<State, Double> cacheC2Values = new TreeMap<>();
		cacheC1Values = findsCS.cacheC1Values;
		cacheC2Values = findsCS.cacheC2Values;
		double simsCSFinalValue = simuation.simulatesCS(initialState, optsCS, cacheC1Values, cacheC2Values, overheadCost, maxOrderQuantity, fixOrderCost, variCost);
		double gap1 = (finalCash -simsCSFinalValue)/finalCash;
		double gap2 = (simFinalValue -simsCSFinalValue)/simFinalValue;	
		System.out.printf("Optimality gap is: %.2f%% or %.2f%%\n", gap1 * 100, gap2 * 100);

		/*******************************************************************
		 * Drawing x Q
		 */
		int minInventorys = 0;
		int maxInventorys = 50; // for drawing pictures
		int minCash = 0;
		int maxCash = 50;
		int xLength = maxCash - minCash + 1;
		int yLength = maxInventorys - minInventorys + 1;
		Drawing drawing = new Drawing();		
		double initialInventory = 3;
		
//		double[][] xQ = new double[xLength][2];
//		int index = 0;
//		for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
//		//for (double initialCash = minCash; initialCash <= maxCash; initialCash++) {
//			period = 1;
//			xQ[index][0] = initialInventory;
//			recursion.getExpectedValue(new CashState(period, initialInventory, iniCash));
//			xQ[index][1] = recursion.getAction(new CashState(period, initialInventory, iniCash));
//			index++;
//		}		
		//drawing.drawXQ(xQ);

		/*******************************************************************
		 * Drawing y G since comupteIfAbsent, we need initializing a new class to draw
		 * Gy; if not, java would not compute sdp again, we must redefine
		 * stateTransition function and immediate Function;
		 */
		// immediate value
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue2 = (state, action, randomDemand) -> {
			double revenue = 0;
			double fixedCost = 0;
			double variableCost = 0;
			double inventoryLevel = 0;
			if (isForDrawGy == true && state.getPeriod() == 1) {
				revenue = price * Math.min(state.getIniInventory(), randomDemand);
				fixedCost = 0;
				variableCost = variCost * state.getIniInventory();
				inventoryLevel = state.getIniInventory() - randomDemand;
			} else {
				revenue = price * Math.min(state.getIniInventory() + action, randomDemand);
				fixedCost = action > 0 ? fixOrderCost : 0;
				variableCost = variCost * action;
				inventoryLevel = state.getIniInventory() + action - randomDemand;
			}
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue - fixedCost - variableCost - holdCosts - overheadCost;
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			return cashIncrement;
		};

		// state transition function 2, for GB
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition2 = (state, action,
				randomDemand) -> {
			double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
					: state.getIniInventory() + action - randomDemand;
			double nextCash = state.getIniCash() + immediateValue2.apply(state, action, randomDemand);
			if (isForDrawGy == true && state.getPeriod() == 1)
				nextCash -= fixOrderCost;
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};
		
		
		CashRecursion recursion2 = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition2,
				immediateValue2, discountFactor);
		double[][] yG2 = new double[xLength][2];
		int index = 0;
		//for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
		for (double initialCash = minCash; initialCash <= maxCash; initialCash++) {
			yG2[index][0] = initialCash; // initialInventory
			yG2[index][1] = recursion2.getExpectedValue(new CashState(period, initialInventory, initialCash)); // iniCash
			index++;
		}
		//CheckKConvexity.check(yG2, fixOrderCost);
		drawing.drawSimpleG(yG2, iniCash, "K transfered in cash GB"); // GB
		
		/*******************************************************************
		 * Drawing another G() that has no fixed ordering cost transition in the 
		 * first period, GA
		 * the difference lies in state transition function
		 */
		
		// state transition function 3
		StateTransitionFunction<CashState, Double, Double, CashState> stateTransition3 = (state, action,
				randomDemand) -> {
			double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
						: state.getIniInventory() + action - randomDemand;
			double nextCash = state.getIniCash() + immediateValue2.apply(state, action, randomDemand);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
		};

		CashRecursion recursion3 = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition3,
				immediateValue2, discountFactor);
		double[][] yG3 = new double[xLength * yLength][3];
		index = 0;
		//for (int initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {
		for (double initialCash = minCash; initialCash <= maxCash; initialCash++) {
			for (initialInventory = minInventoryState; initialInventory <= maxInventorys; initialInventory++) {
			yG3[index][0] = initialCash; // initialInventory
			yG3[index][1] = initialInventory; // initialInventory
			yG3[index][2] = recursion3.getExpectedValue(new CashState(period, initialInventory, initialCash)); // iniCash
			index++;
			}
		}
		WriteToExcel wr = new WriteToExcel();
		wr.writeArrayToExcel(yG3, "GA.xls");
		drawing.drawSimpleG(yG3, iniCash, "K not transfered in cash GA");
		drawing.drawTwoGR(yG3, yG2, initialInventory); //drawing.drawTwoG(yG3, yG2, iniCash);
		
		double[] interPoint = drawing.intersectionPoint(yG3, yG2, iniCash);
		String fileName= "interSectionPoints.xls";
		wr.writeToExcelAppend(interPoint, fileName);
	}
}
