/**
 * @date: Jun 4, 2020
 */
package cash.strongconstraint;

import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import cash.strongconstraint.FindsCS.FindCCrieria;
import sdp.cash.CashRecursion;
import sdp.cash.CashState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 4, 2020
 * @Desc: check the monotony property of GB(y*) - GA(x)- K , GA(x)
 *
 */
public class CheakMonotony {
	public static void main(String[] args) {
		
		double[][] iniMeanDemands = { { 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 },
				{ 21.15, 18.9, 17.7, 16.5, 15.15, 13.95, 12.75, 11.55, 10.35, 9.15}, { 6.6, 9.3, 11.1, 12.9, 16.8, 21.6, 24, 26.4, 31.5, 33.9 },
				{ 12.1,10,7.9,7,7.9,10,12.1,13,12.1,10 }, {15.7,10,4.3,2,4.3,10,15.7,18,15.7,10 },
				{ 41.8, 6.6, 2, 21.8, 44.8, 9.6, 2.6, 17, 30, 35.4 }, { 4.08, 12.16, 37.36, 21.44, 39.12, 35.68, 19.84, 22.48, 29.04, 12.4 },
				{ 4.7, 8.1, 23.6, 39.4, 16.4, 28.7, 50.8, 39.1, 75.4, 69.4 }, { 4.4, 11.6, 26.4, 14.4, 14.6, 19.8, 7.4, 18.3, 20.4, 11.4 },
				{ 4.9, 18.8, 6.4, 27.9, 45.3, 22.4, 22.3, 51.7, 29.1, 54.7 } };
		
		
		double[] K = {10, 15, 20};
		double[] v = {1};
		double[] B0 = { 3, 5, 7}; // ini cash can order 4 or 6 items
		double[] p = { 5, 6, 7};  // margin is 4, 5, 6
		double[] h = {0};
		double salvageValue = 0;	
		
		FindCCrieria criteria = FindCCrieria.XRELATE;
		double truncationQuantile = 0.999;
		int stepSize = 1;
		double overheadCost = 0; // minimum cash balance the retailer can withstand
		double maxOrderQuantity = 150; // maximum ordering quantity when having enough cash
		double minInventoryState = 0;
		double maxInventoryState = 200;
		double minCashState = -100; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 1500;
		double discountFactor = 1; // generally 1 for the cash constrained problem
		boolean isForDrawGy = true;
		
		/*******************************************************************
		 * set demands length, for testing 
		 */
		int newLength = 3;
		double[][] meanDemands = new double[iniMeanDemands.length][newLength];
		for (int i = 0; i < iniMeanDemands.length; i++)
			for (int j = 0; j < newLength; j++) {
				meanDemands[i][j] = iniMeanDemands[i][j];
			}
		
		
		for (int idemand = 0; idemand < meanDemands.length; idemand++)
			for (int iK = 0; iK < K.length; iK++)
				for (int iv = 0; iv < v.length; iv++)
					for (int ip = 0; ip < p.length; ip++)
						for (int ih = 0; ih < h.length; ih++)
							for (int iB = 0; iB < B0.length; iB++) {
								double[] meanDemand = meanDemands[idemand];
								double fixOrderCost = K[iK];
								double variCost = v[iv];
								double price = p[ip];
								double iniCash = fixOrderCost + overheadCost + variCost * B0[iB];
								double holdingCost = h[ih];
								double iniInventory = 0;
								
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
								int maxInventorys = 30; // for drawing pictures
								int minCash = (int) fixOrderCost;
								int maxCash = (int) fixOrderCost + 30;
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
								
								for (double initialCash = minCash; initialCash <= maxCash; initialCash++) {
									columnIndex = 0;
									for (double initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {			
										yG[index][0] = initialCash; // initialInventory
										yG[index][1] = initialInventory;
										yG[index][2] = recursion.getExpectedValue(new CashState(period, initialInventory, initialCash)); // iniCash
										resultTableF[rowIndex][columnIndex] = yG[index][2];
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

//								// state transition function, for GB
//								StateTransitionFunction<CashState, Double, Double, CashState> stateTransition2 = (state, action,
//										randomDemand) -> {
//									double nextInventory = isForDrawGy && state.getPeriod() == 1 ? state.getIniInventory() - randomDemand
//											: state.getIniInventory() + action - randomDemand;
//									double nextCash = state.getIniCash() + immediateValue2.apply(state, action, randomDemand);
//									if (isForDrawGy == true && state.getPeriod() == 1)
//										nextCash -= fixOrderCost;
//									nextCash = nextCash > maxCashState ? maxCashState : nextCash;
//									nextCash = nextCash < minCashState ? minCashState : nextCash;
//									nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
//									nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
//									return new CashState(state.getPeriod() + 1, nextInventory, nextCash);
//								};
								
								

								
//								CashRecursion recursion2 = new CashRecursion(OptDirection.MAX, pmf, getFeasibleAction, stateTransition2,
//										immediateValue2, discountFactor);
//								double[][] yG2 = new double[xLength * RLength][3];
//								double[][] resultTableGB = new double[RLength][xLength];
//								
//								for (double initialCash = minCash; initialCash <= maxCash; initialCash++) {
//									columnIndex = 0;
//									for (double initialInventory = minInventorys; initialInventory <= maxInventorys; initialInventory++) {			
//									yG2[index][0] = initialCash; // initialInventory
//									yG2[index][1] = initialInventory;
//									yG2[index][2] = recursion2.getExpectedValue(new CashState(period, initialInventory, initialCash)); // iniCash
//									resultTableGB[rowIndex][columnIndex] = yG2[index][2];
//									index++;
//									columnIndex++;
//									}
//									rowIndex++;
//								}
								
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
								double[][] resultTableGA = new double[RLength][xLength];
								rowIndex = 0;
								for (double initialCash = minCash; initialCash <= maxCash; initialCash++) {
									columnIndex = 0;
									for (double initialInventory = minInventoryState; initialInventory <= maxInventorys; initialInventory++) {
									yG3[index][0] = initialCash; // initialInventory
									yG3[index][1] = initialInventory; // initialInventory
									yG3[index][2] = recursion3.getExpectedValue(new CashState(period, initialInventory, initialCash)); // iniCash
									resultTableGA[rowIndex][columnIndex] = yG3[index][2];
									index++;
									columnIndex++;
									}
									rowIndex++;
								}
								
								double[][] resultMinusFG = recursion3.getMinusFG(resultTableGA, resultTableF, minCash, fixOrderCost, variCost);
								double time = (System.currentTimeMillis() - currTime) / 1000;
								System.out.println("running time is " + time + "s");

							}
		

	}

}
