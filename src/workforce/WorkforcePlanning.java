package workforce;

import java.util.function.Function;

import sdp.inventory.CheckKConvexity;
import sdp.inventory.Drawing;
import sdp.inventory.Recursion;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.Recursion.OptDirection;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;




/**
 * @author chen
 * @description: optimal ordering quantity for a single period problem is 
 * F^{-1}((\pi-h-v)/\pi)-x+w.
 *
 */
public class WorkforcePlanning {

	public static void main(String[] args) {
		double[] dimissionRate = {0.5, 0.2, 0.5, 0.3};
		int T = dimissionRate.length;
		int iniStaffNum = 80;
		double fixCost = 10000;
		double unitVariCost = 1000;
		double salary = 3000;
		double unitPenalty = 5000;		
		int[] minStaffNum = {90, 90, 90, 90};	
		
		int maxHireNum = 300;
		int stepSize = 1;
		double truncQuantile = 0.999;
		boolean isForDrawGy = true;
		
		int minX = 10;
		int maxX = 300; // for drawing pictures
		
		// feasible actions
		Function<StaffState, int[]> getFeasibleAction = s -> {
			int[] feasibleActions = new int[(maxHireNum / stepSize) + 1];
			int index = 0;
			for (int i = 0; i <= maxHireNum; i = i + stepSize) {
				feasibleActions[index] = i;
				index++;
			}
			return feasibleActions;
		};
		
		// state transition function
		StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition = (state, action, randomDemand) -> {
			int nextStaffNum = state.iniStaffNum + action - randomDemand;
			return new StaffState(state.period + 1, nextStaffNum);
		};
		
		// immediate value function
		ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue = (state, action, randomDemand) -> {
			double fixHireCost = action > 0 ? fixCost: 0;
			double variHireCost = unitVariCost * action;
			int nextStaffNum = state.iniStaffNum + action - randomDemand;
			double salaryCost = salary * nextStaffNum;
			int t = state.period - 1;
			double penaltyCost = nextStaffNum > minStaffNum[t] ? 0 : unitPenalty * (minStaffNum[t] - nextStaffNum);
			double totalCosts = fixHireCost + variHireCost + salaryCost + penaltyCost;			
			return totalCosts;
		};
		
		/*******************************************************************
		 * Solve
		 */
		StaffRecursion recursion = new StaffRecursion(getFeasibleAction, stateTransition, immediateValue, dimissionRate, truncQuantile);
		int period = 1;
		StaffState initialState = new StaffState(period, iniStaffNum);
		long currTime = System.currentTimeMillis();
		double opt = recursion.getExpectedValue(initialState);
		System.out.println("final optimal expected cost is: " + opt);
		System.out.println("optimal hiring number in the first priod is : " + recursion.getAction(initialState));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		double[][] optTable = recursion.getOptTable();
		
		
		/*******************************************************************
		 * Drawing
		 */
		int xLength = maxX - minX + 1;
		double[][] xQ = new double[xLength][2];
		int index = 0;
		for (int initialInventory = minX; initialInventory <= maxX; initialInventory++) {
			period = 1;
			xQ[index][0] = initialInventory;
			recursion.getExpectedValue(new StaffState(period, initialInventory));
			xQ[index][1] = recursion.getAction(new StaffState(period, initialInventory));
			index++;
		}
		Drawing.drawXQ(xQ);
		
		// since comupteIfAbsent, we need initializing a new class to draw Gy; if not, java would not compute sdp again
		// must redefine stateTransition function and immediate Function;
		StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition2 = (state, action, randomDemand) -> {
			int nextStaffNum = isForDrawGy && state.period == 1 ? state.iniStaffNum - randomDemand
					: state.iniStaffNum + action - randomDemand;
			return new StaffState(state.period + 1, nextStaffNum);
		};

		ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue2 = (state, action, randomDemand) -> {
			double fixHireCost = 0;
			double variHireCost = unitVariCost * action;
			int nextStaffNum = state.iniStaffNum + action - randomDemand;
			double salaryCost = salary * nextStaffNum;
			int t = state.period - 1;
			double penaltyCost = nextStaffNum > minStaffNum[t] ? 0 : unitPenalty * (minStaffNum[t] - nextStaffNum);
			double totalCosts = fixHireCost + variHireCost + salaryCost + penaltyCost;			
			return totalCosts;
		};

		StaffRecursion recursion2 = new StaffRecursion(getFeasibleAction, stateTransition2, immediateValue2, dimissionRate, truncQuantile);
		
		double[][] yG = new double[xLength][2];
		index = 0;
		for (int initialStaff = minX; initialStaff <= maxX; initialStaff++) {
			yG[index][0] = initialStaff;
			yG[index][1] = recursion2.getExpectedValue(new StaffState(period, initialStaff));
			index++;
		}
		CheckKConvexity.check(yG, fixCost);
		Drawing.drawSimpleG(yG);
		Drawing.drawGAndsS(yG, fixCost);
		
//		StaffRecursion recursion2 = new StaffRecursion(getFeasibleAction, stateTransition, immediateValue, dimissionRate, truncQuantile);
//		double opt2 = recursion2.getExpectedValueNoHireFirst(initialState);
//		System.out.println("final optimal expected cost if not hiring in the first period is: " + opt2);
	}

}
