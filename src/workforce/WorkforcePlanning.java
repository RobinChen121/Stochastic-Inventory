package workforce;

import java.util.Arrays;
import java.util.function.Function;

import milp.MIPWorkforce;
import sdp.inventory.CheckKConvexity;
import sdp.inventory.Drawing;
import sdp.inventory.FitsS;
import sdp.inventory.Recursion;
import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.Recursion.OptDirection;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.BinomialDist;
import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;




/**
 * @author chen
 * @description: optimal ordering quantity for a single period problem is 
 * F^{-1}(((\pi-h)(1-p)-v)/(\pi(1-p)))+w = y*.
 *
 *not only (s, S) policy optimal, (R, S) may also be optimal 
 */
public class WorkforcePlanning {

	public static void main(String[] args) {
		double[] turnoverRate;
		turnoverRate = new double[] {0.3, 0.9, 0.3};
		//Arrays.fill(turnoverRate, 0.5);
		int T = turnoverRate.length;
		
		int iniStaffNum = 0;
		double fixCost = 100;
		double unitVariCost = 0;
		double salary = 70;
		double unitPenalty = 80;		
		int[] minStaffNum = {30, 50, 40};	
		
		int maxHireNum = 500;
		int maxX = maxHireNum; // for drawing pictures
		int stepSize = 1;
		boolean isForDrawGy = true;
		int segmentNum = 10;
		
		int minX = 0;
		int xLength = maxX - minX + 1;
		
		// get pmf for every possible x
		double[][][][] pmf = new double[T][][][];
		for (int t = 0; t < T; t++) {
			pmf[t] = new double[xLength][][];
			for (int i = minX; i <= maxX; i++) {
				pmf[t][i] = new double[i+1][2];
				if (i == 0) {
					pmf[t][0] = new double[][]{{0, 1}};
				}
				else {			
					BinomialDist distribution = new BinomialDist(i, turnoverRate[t]);
					for (int j = 0; j < i+1; j++) {
						pmf[t][i][j][0] = j;		
						pmf[t][i][j][1] = distribution.prob(j);
					}
				}			
			}		
		}
		
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
			nextStaffNum  = nextStaffNum  > maxX ? maxX : nextStaffNum;
			nextStaffNum  = nextStaffNum  < minX ? minX : nextStaffNum;
			return new StaffState(state.period + 1, nextStaffNum);
		};
		
		// immediate value function
		ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue = (state, action, randomDemand) -> {
			double fixHireCost = action > 0 ? fixCost : 0;
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
		StaffRecursion recursion = new StaffRecursion(getFeasibleAction, stateTransition, immediateValue, pmf, T);
		int period = 1;
		StaffState initialState = new StaffState(period, iniStaffNum);
		long currTime = System.currentTimeMillis();
		double opt = recursion.getExpectedValue(initialState);
		System.out.println("final optimal expected cost is: " + opt);
		int optQ = recursion.getAction(initialState);
		System.out.println("optimal hiring number in the first priod is : " + optQ);
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		double[][] optTable = recursion.getOptTable();
		
		
		/*******************************************************************
		 * find s and S from SDP and simulate.
		 * when >= s, not order.
		 */
		FitsS findsS = new FitsS(Integer.MAX_VALUE, T);
		double[][] optsS = findsS.getSinglesS(optTable);
		System.out.println("single s, S level: " + Arrays.deepToString(optsS));
		
		SimulatesS simulate = new SimulatesS(recursion, T, turnoverRate);
		double sim  = simulate.simulatesS(initialState, optsS);
		
		System.out.printf("simulated value is %.2f\n", sim);
		System.out.printf("simulated gap is %.2f%%\n", (sim - opt)*100/opt);
		
		
		/*******************************************************************
		 * piecewise MIP
		 */
		System.out.println("**********************************************");
		MIPWorkforce mip = new MIPWorkforce(iniStaffNum, fixCost, unitVariCost, salary, unitPenalty, minStaffNum, turnoverRate);
		
		double mipObj = mip.pieceApprox(segmentNum);
		double[][] sS = mip.getsS(segmentNum);
		System.out.println("s, S by mip are: " + Arrays.deepToString(sS));
		double sim2 = simulate.simulatesS(initialState, sS);

		System.out.printf("simulated value is %.2f\n", sim2);
		System.out.printf("simulated gap is %.2f%%\n", (sim2 - opt)*100/opt);
		
		
		/*******************************************************************
		 * Drawing
		 * xQ
		 */
//		double[][] xQ = new double[xLength][2];
//		int index = 0;
//		for (int initialInventory = minX; initialInventory <= maxX; initialInventory++) {
//			period = 1;
//			xQ[index][0] = initialInventory;
//			recursion.getExpectedValue(new StaffState(period, initialInventory));
//			xQ[index][1] = recursion.getAction(new StaffState(period, initialInventory));
//			index++;
//		}
//		Drawing.drawXQ(xQ);
//		
//		/*******************************************************************
//		 * Drawing
//		 * G(y, x), G should also related with x since x affected the demand distribution.
//		 * since comupteIfAbsent, we need initializing a new class to draw Gy; if not, java would not compute sdp again.
//		 * must redefine stateTransition function and immediate Function.
//		 */
//		StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition2 = (state, action, randomDemand) -> {
//			int nextStaffNum = state.period == 1 ? state.iniStaffNum - randomDemand : state.iniStaffNum + action - randomDemand;
//			return new StaffState(state.period + 1, nextStaffNum);
//		};
//
//		ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue2 = (state, action, randomDemand) -> {
//			double fixHireCost;
//			double variHireCost;
//			int nextStaffNum;
//			if (state.period == 1) {
//				fixHireCost = 0;
//				variHireCost = unitVariCost * state.iniStaffNum;
//				nextStaffNum = state.iniStaffNum - randomDemand;
//			}
//			else {
//				fixHireCost = action > 0 ? fixCost : 0;
//				variHireCost = unitVariCost * action;
//				nextStaffNum = state.iniStaffNum + action - randomDemand;
//			}
//			double salaryCost = salary * nextStaffNum;
//			int t = state.period - 1;
//			double penaltyCost = nextStaffNum > minStaffNum[t] ? 0 : unitPenalty * (minStaffNum[t] - nextStaffNum);
//			double totalCosts = fixHireCost + variHireCost + salaryCost + penaltyCost;			
//			return totalCosts;
//		};
//
//		StaffRecursion recursion2 = new StaffRecursion(getFeasibleAction, stateTransition2, immediateValue2, pmf, T);
//		
//		double[][] yG = new double[xLength][2];
//		index = 0;
//		for (int initialStaff = minX; initialStaff <= maxX; initialStaff++) {
//			yG[index][0] = initialStaff;
//			yG[index][1] = recursion2.getExpectedValue(new StaffState(period, initialStaff), iniStaffNum);
//			index++;
//		}
//		CheckKConvexity CheckK = new CheckKConvexity();
//		CheckK.check(yG, fixCost);
//		Drawing.drawSimpleG(yG);
//		Drawing.drawGAndsS(yG, fixCost);
		
//		StaffRecursion recursion2 = new StaffRecursion(getFeasibleAction, stateTransition, immediateValue, turnoverRate, truncQuantile);
//		double opt2 = recursion2.getExpectedValueNoHireFirst(initialState);
//		System.out.println("final optimal expected cost if not hiring in the first period is: " + opt2);
		
//		BinomialDist dist = new BinomialDist(optQ, turnoverRate[0]);
//		int QStar = dist.inverseFInt((unitPenalty-salary*(1-turnoverRate[0])-unitVariCost)/unitPenalty) - iniStaffNum + minStaffNum[0];
//		System.out.println("the optimal ordering quantity in the first period is: " + QStar);
	}

}
