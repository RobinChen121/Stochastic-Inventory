package workforce;

import java.util.concurrent.CompletionStage;
import java.util.function.Function;

import sdp.inventory.Recursion;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.Recursion.OptDirection;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;

public class WorkforcePlanning {

	public static void main(String[] args) {
		double[] dimissionRate = {0.1};
		int T = dimissionRate.length;
		int iniStaffNum = 100;
		double fixCost = 1000;
		double unitVariCost = 100;
		double salary = 5000;
		double unitPenalty = 50;
		
		int minStaffNum = 90;	
		int maxHireNum = 100;
		int stepSize = 1;
		double truncQuantile = 0.99;
		
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
			double penaltyCost = nextStaffNum > minStaffNum ? 0 : unitPenalty * (minStaffNum - nextStaffNum);
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
	}

}
