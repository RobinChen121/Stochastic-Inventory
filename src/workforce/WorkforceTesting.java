package workforce;

import java.util.Arrays;
import java.util.function.Function;

import milp.MIPWorkforce;
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
		double[] salarys = {30, 100, 200};
		double[] unitPenaltys = {50, 250, 2200};
		int[][] minStaffs = {{40, 40, 40, 40, 40, 40, 40, 40},
				{10, 20, 30, 40, 50, 60, 70, 80},
				{80, 70, 60, 50, 40, 30, 20, 10},
				{40, 50, 80, 60, 40, 50, 80, 60},
		};
		int segmentNum = 14;
		
		WriteToExcelTxt wr = new WriteToExcelTxt();
		String fileName = "results.xls";
		String headString =  
				"turnoverRate" + "\t" + "fixCost" + "\t"  + "salary" + "\t" + "unitPenalty" + "\t" + "iMinStaff" + "\t" +
		         "Q*" + "\t" + "ExpectedCosts"  + "\t" + "time"+ "\t" + "simCosts"+ "\t" + "gapPercent" +"\t" +
						    "MIPcost" + "\t" + "timeMip" + "\t" + "gapMIP" + "\t" + "MIPsScost" + "\t" + "timesS" + "\t" + "gapMIPsS";
		wr.writeToFile(fileName, headString);
		
		int m = turnoverRates.length;
		for (int iRate = 1; iRate < 2; iRate++)
			for (int iFix = 0; iFix < m; iFix++)
				for (int iSalary = 0; iSalary < m; iSalary ++)
					for (int iPenalty = 2; iPenalty < 3; iPenalty ++) 
						for (int iMinStaff = 0; iMinStaff < m + 1; iMinStaff ++) {
			
			int T = 8;
			double[] turnoverRate = new double[T];
			Arrays.fill(turnoverRate, turnoverRates[iRate]);
			int iniStaffNum = 0;
			double fixCost = fixCosts[iFix];
			double unitVariCost = 0;
			double salary = salarys[iSalary];
			double unitPenalty = unitPenaltys[iPenalty];		
			int[] minStaffNum = minStaffs[iMinStaff];	
			
			int maxHireNum = 800;
			int stepSize = 1;
			boolean isForDrawGy = true;
			
			int minX = iniStaffNum;
			int maxX = maxHireNum; // for drawing pictures
			
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
			
			/*******************************************************************
			 * find s and S from MIP and simulate.
			 * when >= s, not order.
			 */
			System.out.println("**********************************************");
			MIPWorkforce mip = new MIPWorkforce(iniStaffNum, fixCost, unitVariCost, salary, unitPenalty, minStaffNum, turnoverRate);
			
			currTime = System.currentTimeMillis();
			double mipObj = mip.pieceApprox(segmentNum);
			double timeMip = (System.currentTimeMillis() - currTime) / 1000;
			System.out.println("running time for mip is " + timeMip + "s");
			System.out.printf("mip value for expected cost is %.2f\n", mipObj);
			double gapMip = (mipObj - opt)*100/opt;
			System.out.printf("gap for mip  is %.2f%%\n", gapMip);
			
			currTime = System.currentTimeMillis();
			double[][] sS = mip.getsS(segmentNum);
			double timeMipsS = (System.currentTimeMillis() - currTime) / 1000;
			System.out.println("running time for mip is " + timeMipsS + "s");
			
			System.out.println("s, S by mip are: " + Arrays.deepToString(sS));
			double sim2 = simulate.simulatesS(initialState, sS);

			System.out.printf("simulated value for mip sS is %.2f\n", sim2);
			double gapsS = (sim2 - opt)*100/opt;
			System.out.printf("simulated gap for mip sS is %.2f%%\n", gapsS);
			System.out.println("**********************************************");
			System.out.println("**********************************************");
			System.out.println("**********************************************");
			
			double[] out = new double[]{turnoverRate[0], fixCost, salary, unitPenalty, iMinStaff, optQ, opt, time, sim, gapPercent,
					mipObj, timeMip, gapMip, sim2, timeMipsS, gapsS};
			wr.writeToExcelAppend(out, fileName);		
		}
	}

}
