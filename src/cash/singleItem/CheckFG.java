/**
 * @date: Jun 7, 2020
 */
package cash.singleItem;

import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.cash.CashRecursion;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToCsv;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 7, 2020
 * @Desc: check the values and relation of F and G
 *
 */
public class CheckFG {

	public static void main(String[] args) {
		double[] meanDemand = {8, 8, 8};
		double iniCash = 13;
		double iniInventory = 0;
		double fixOrderCost = 10;
		double variCost = 1;
		double price = 8;
		double salvageValue = 0.5;
		double holdingCost = 0;		
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
				//.mapToObj(i -> new NormalDist(meanDemand[i], 0.25 * meanDemand[i])) // can be changed to other distributions
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
		 * Solve F(x, R)
		 */
		int minInventorys = 0;
		int maxInventorys = 100; // for drawing pictures
		int minCash = 0;
		int maxCash = (int) fixOrderCost + 80;
		int RLength = maxCash - minCash + 1;
		int xLength = maxInventorys - minInventorys + 1;
		int period = 1;
		int index = 0;
		int rowIndex = 0;
		int columnIndex = 0;
		long currTime = System.currentTimeMillis();
		
		CashRecursion recursion = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition,
				immediateValue, discountFactor);
		
		double[][] yG = new double[xLength * RLength][3];
		double[][] resultTableF = new double[RLength][xLength];
		double[][] resultTableQ = new double[RLength][xLength];
		
		for (double initialCash = minCash; initialCash <= maxCash; initialCash++) {
			columnIndex = 0;
			for (double initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {			
				yG[index][0] = initialCash; // initialInventory
				yG[index][1] = initialInventory;
				yG[index][2] = recursion.getExpectedValue(new CashState(period, initialInventory, initialCash)); // iniCash
				resultTableF[rowIndex][columnIndex] = yG[index][2];
				resultTableQ[rowIndex][columnIndex] = recursion.getAction(new CashState(period, initialInventory, initialCash));
				index++;
				columnIndex++;
			}
			rowIndex++;
		}				
		
		// immediate value for GB and GA
		ImmediateValueFunction<CashState, Double, Double, Double> immediateValue2 = (state, action, randomDemand) -> {
			double revenue = 0;
			double fixedCost = 0;
			double variableCost = 0;
			double inventoryLevel = 0;
			if (isForDrawGy == true && state.getPeriod() == 1) {
				revenue = price * Math.min(state.getIniInventory(), randomDemand);
				fixedCost = 0;
				variableCost = 0; // variCost * state.getIniInventory(); // be careful
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
		
		// state transition function for GA
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
		double[][] yG3 = new double[xLength * RLength][3];
		index = 0;
		double[][] resultTableGA = new double[xLength][RLength];
		rowIndex = 0;
		for (double initialInventory = minInventoryState; initialInventory <= maxInventorys; initialInventory++) {
			columnIndex = 0;
			for (double initialCash = minCash; initialCash <= maxCash; initialCash++) {
				yG3[index][0] = initialInventory; // initialInventory
				yG3[index][1] = initialCash; // initialCash
				yG3[index][2] = recursion3.getExpectedValue(new CashState(period, initialInventory, initialCash)) - variCost * initialInventory; 
				resultTableGA[rowIndex][columnIndex] = yG3[index][2]; // careful, minus cy
				index++;
				columnIndex++;
			}
			rowIndex++;
		}
		
		double[][] resultH= recursion3.getH(resultTableGA, minCash, fixOrderCost, variCost);
		
		System.out.println("Single-crossing property: " + recursion3.checkSingleCrossing(resultH));
		
		WriteToCsv wr = new WriteToCsv();
		wr.writeArrayCSVLabel(resultTableQ, minCash, minInventorys, "Q.csv");
		wr.writeArrayCSVLabel(resultTableF, minCash, minInventorys, "F.csv");
		wr.writeArrayCSVLabel(resultTableGA, minCash, minInventorys,  "G.csv");
		wr.writeArrayCSVLabel(resultH, minCash, minInventorys,  "H.csv");
//		wr.writeArrayCSV(recursion3.getH3Column(resultTableGA, minCash, fixOrderCost, variCost), "G3column.csv");
		
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
	}

}
