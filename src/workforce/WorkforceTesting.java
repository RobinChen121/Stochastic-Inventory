package workforce;

import java.util.Arrays;
import java.util.function.Function;

import sdp.inventory.CheckKConvexity;
import sdp.inventory.Drawing;
import sdp.inventory.FitsS;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.probdist.BinomialDist;

public class WorkforceTesting {
	public static void main(String[] args) {
		double[] turnoverRates = {0.1, 0.5, 0.9};
		double[] fixCosts = {50, 1000, 2000};
		double[] salarys = {10, 100, 200};
		double[] unitPenaltys = {50, 250, 2200};
		
		WriteToExcelTxt wr = new WriteToExcelTxt();
		String fileName = "results.xls";
		String headString =  
				"turnoverRate" + "\t" + "fixCost" + "\t"  + "salary" + "\t" + "unitPenalty" + "\t" +
		         "Q*" + "\t" + "ExpectedCosts"  + "\t" + "time"+ "\t" + "simCosts"+ "\t" + "gapPercent";
		wr.writeToFile(fileName, headString);
		
		int m = turnoverRates.length;
		for (int iRate = 1; iRate < 2; iRate++)
			for (int iFix = 0; iFix < 1; iFix++)
				for (int iSalary = 0; iSalary < 1; iSalary ++)
					for (int iPenalty = 2; iPenalty < m; iPenalty ++) {
			
			int T = 12;
			double[] turnoverRate = new double[T];
			Arrays.fill(turnoverRate, turnoverRates[iRate]);
			int iniStaffNum = 0;
			double fixCost = fixCosts[iFix];
			double unitVariCost = 0;
			double salary = salarys[iSalary];
			double unitPenalty = unitPenaltys[iPenalty];		
			int[] minStaffNum = {80, 10, 58, 65, 15, 30, 45, 60, 10, 56, 49, 70};	
			
			int maxHireNum = 300;
			int stepSize = 1;
			boolean isForDrawGy = true;
			
			int minX = iniStaffNum;
			int maxX = 300; // for drawing pictures
			
			// get pmf for every possible x
			int xLength = maxX - minX + 1;
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
			double optQ = recursion.getAction(initialState);
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
			double gapPercent = (sim - opt)*100/opt;
			System.out.printf("simulated gap is %.2f%%\n", gapPercent);
			System.out.println("****************************************************");
			System.out.println();
						
			/*******************************************************************
			 * Drawing
			 * G(y, x), G should also related with x since x affected the demand distribution.
			 * since comupteIfAbsent, we need initializing a new class to draw Gy; if not, java would not compute sdp again.
			 * must redefine stateTransition function and immediate Function.
			 */
//			StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition2 = (state, action, randomDemand) -> {
//				int nextStaffNum = state.period == 1 ? state.iniStaffNum - randomDemand : state.iniStaffNum + action - randomDemand;
//				return new StaffState(state.period + 1, nextStaffNum);
//			};
//	
//			ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue2 = (state, action, randomDemand) -> {
//				double fixHireCost;
//				double variHireCost;
//				int nextStaffNum;
//				if (state.period == 1) {
//					fixHireCost = 0;
//					variHireCost = unitVariCost * state.iniStaffNum;
//					nextStaffNum = state.iniStaffNum - randomDemand;
//				}
//				else {
//					fixHireCost = action > 0 ? fixCost : 0;
//					variHireCost = unitVariCost * action;
//					nextStaffNum = state.iniStaffNum + action - randomDemand;
//				}
//				double salaryCost = salary * nextStaffNum;
//				int t = state.period - 1;
//				double penaltyCost = nextStaffNum > minStaffNum[t] ? 0 : unitPenalty * (minStaffNum[t] - nextStaffNum);
//				double totalCosts = fixHireCost + variHireCost + salaryCost + penaltyCost;			
//				return totalCosts;
//			};
	
//			StaffRecursion recursion2 = new StaffRecursion(getFeasibleAction, stateTransition2, immediateValue2, pmf, T);
//			
//			
//			double[][] yG = new double[xLength][2];
//			int index = 0;
//			for (int initialStaff = minX; initialStaff <= maxX; initialStaff++) {
//				yG[index][0] = initialStaff;
//				yG[index][1] = recursion2.getExpectedValue(new StaffState(period, initialStaff), iniStaffNum);
//				index++;
//			}
//			CheckKConvexity CheckK = new CheckKConvexity();
//			Boolean ck = true; //CheckK.check(yG, fixCost);
			
			double[] out = new double[]{turnoverRate[0], fixCost, salary, unitPenalty, optQ, opt, time, sim, gapPercent};
			wr.writeToExcelAppend(out, fileName);
			

		}
	}

}
