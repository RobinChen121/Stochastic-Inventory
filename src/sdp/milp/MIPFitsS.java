package sdp.milp;

import java.util.ArrayList;
import java.util.Arrays;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 12, 2018---8:48:27 PM
 * @description: find s, S for capacitated lot sizing problem by fitting via Cplex
 * 
 * @note: in this class, you need have a cplex.jar library, download CPLEX software from 
 * https://www.ibm.com/analytics/cplex-optimizer
 */

public class MIPFitsS {
	
	int maxOrderQuantity;
	int T;
	
	public MIPFitsS(int maxOrderQuantity, int T) {
		this.maxOrderQuantity = maxOrderQuantity;
		this.T = T;
	}
	

	/**
	 * 
	 * @param optTable   [t, i, Q] in each row
	 * @param maxOrderQuantity  maximum ordreing quantity
	 * @return return s level index in one period
	 */
	public int[] levelIndex(double[][] optTable) {
		ArrayList<Integer> indexArr = new ArrayList<>();
		boolean mark = false;
		for (int j = 0; j < optTable.length; j++) {
			if (optTable[j][2] < maxOrderQuantity && !mark) {
				mark = true;
			} 
			else if (optTable[j][2] == maxOrderQuantity && mark && j != optTable.length - 1) {
				mark = false;
				indexArr.add(j);
			}
			if (optTable[j][2] == 0) {
				indexArr.add(j);
				break;
			}
			if (j == optTable.length - 1) {
				indexArr.add(j);
			}
		}
		return indexArr.stream().mapToInt(p -> p.intValue()).toArray();
	}

	/**
	 * 
	 * @param lb lower value bound
	 * @param upIndex up value bound
	 * @param tOptTable optimal sdp table in a period, 
	 *        [t, i, Q] in each row
	 * @return a fitting value for some
	 */
	public double minSquare(double lb, int upIndex, double[][] tOptTable) {
		int lowIndex = 0;
		double realS;
		for (int i = 0; i < tOptTable.length; i++)
			if (tOptTable[i][2] != maxOrderQuantity) {
				lowIndex = i;
				break;
			}
		try {
			IloCplex cplex = new IloCplex();
			double ub = 10000;
			IloNumVar x = cplex.numVar(lb, ub);
			realS = tOptTable[lowIndex][1] + tOptTable[lowIndex][2];
			IloNumExpr obj = cplex.prod(cplex.diff(x, realS), cplex.diff(x, realS));
			for (int i = lowIndex + 1; i <= upIndex; i++) {
				if (tOptTable[i][2] != maxOrderQuantity) {
					realS = tOptTable[i][1] + tOptTable[i][2];
					obj = cplex.sum(obj, cplex.prod(cplex.diff(x, realS), cplex.diff(x, realS)));
				}
			}
			cplex.addMinimize(obj);
			if (cplex.solve()) {
				return cplex.getValue(x);
			}
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		return 0;
	}
	
	public double[][] getSinglesS(double[][] optimalTable) {
		double[][] optimalsS = new double[T][2];
		optimalsS[0][0] = optimalTable[0][1] + 1;
		optimalsS[0][1] = optimalTable[0][1] + optimalTable[0][2];
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i)
					.map(p -> Arrays.stream(p).toArray()).toArray(double[][]::new);
			int[] numIndex = levelIndex(tOptTable);
			if (numIndex.length == 1 && numIndex[0] != 0) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				if (numIndex[0] == tOptTable.length - 1 && tOptTable[numIndex[0]][2] == maxOrderQuantity) {
					optimalsS[t][0] = tOptTable[numIndex[0]][1] + 1;
					optimalsS[t][1] = tOptTable[numIndex[0]][1] + tOptTable[numIndex[0]][2];
				}
				
			} else if (numIndex.length == 1 && numIndex[0] == 0) {  // s, S are both zeros
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0]][1];
			} else if (numIndex.length == 0) { // order at max
				optimalsS[t][0] = tOptTable[tOptTable.length - 1][1];
				optimalsS[t][1] = maxOrderQuantity * 10;
			} else {
				int SIndex = numIndex[numIndex.length - 1]; // 最后一个拐点之前的全部 S 数据进行拟合
				optimalsS[t][0] = tOptTable[SIndex][1];
				optimalsS[t][1] = minSquare(optimalsS[t][0], SIndex, tOptTable);
			}
		}
		return optimalsS;
	}
	
	public double[][] getTwosS(double[][] optimalTable){
		double[][] optimalsS = new double[T][4];
		optimalsS[0][0] = optimalTable[0][1] + 1;
		optimalsS[0][1] = optimalTable[0][1] + optimalTable[0][2];
		optimalsS[0][2] = optimalsS[0][0]; 
		optimalsS[0][3] = optimalTable[0][1] + optimalTable[0][2];
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
											.toArray(double[][] :: new);
			int[] numIndex = levelIndex(tOptTable);
			if (numIndex.length == 1 && numIndex[0] != 0) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				if (numIndex[0] == tOptTable.length - 1 && tOptTable[numIndex[0]][2] == maxOrderQuantity) {
					optimalsS[t][0] = tOptTable[numIndex[0]][1] + 1;
					optimalsS[t][1] = tOptTable[numIndex[0]][1] + tOptTable[numIndex[0]][2];
				}
				optimalsS[t][2] = optimalsS[t][0];
				optimalsS[t][3] = optimalsS[t][1];
				
			}
			else if (numIndex.length == 1 && numIndex[0] == 0) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0]][1];
				optimalsS[t][2] = tOptTable[numIndex[0]][1];
				optimalsS[t][3] = tOptTable[numIndex[0]][1];
			}
			else if (numIndex.length == 0) {
				optimalsS[t][0] = tOptTable[tOptTable.length - 1][1];
				optimalsS[t][1] = maxOrderQuantity * 10;
				optimalsS[t][2] = tOptTable[tOptTable.length - 1][1];
				optimalsS[t][3] = maxOrderQuantity * 10;
			}
			else if (numIndex.length == 2) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				optimalsS[t][2] = tOptTable[numIndex[1]][1];
				optimalsS[t][3] = tOptTable[numIndex[1] - 1][1] + tOptTable[numIndex[1] - 1][2];
				if (numIndex[1] == tOptTable.length - 1 && tOptTable[numIndex[1]][2] == maxOrderQuantity) {
					optimalsS[t][2] = tOptTable[numIndex[1]][1] + 1;
					optimalsS[t][3] = tOptTable[numIndex[1]][1] + tOptTable[numIndex[1]][2];
				}
			}
			
			else {
				int sIndex2 = numIndex[numIndex.length - 1];
				optimalsS[t][2] = tOptTable[sIndex2][1];
				optimalsS[t][3] = tOptTable[sIndex2 - 1][1] + tOptTable[sIndex2 - 1][2];
				int sIndex1 = numIndex[numIndex.length - 2];
				optimalsS[t][0] = tOptTable[sIndex1][1];
				optimalsS[t][1] = minSquare(optimalsS[t][0], sIndex1, tOptTable);
			}
		}
		//System.out.println(Arrays.deepToString(optimalsS));
		return optimalsS;
	}
	
	public double[][] getThreesS(double[][] optimalTable){
		double[][] optimalsS = new double[T][6];
		optimalsS[0][0] = optimalTable[0][1] + 1;
		optimalsS[0][1] = optimalTable[0][1] + optimalTable[0][2];
		optimalsS[0][2] = optimalsS[0][0];
		optimalsS[0][3] = optimalTable[0][1] + optimalTable[0][2];
		optimalsS[0][4] = optimalsS[0][0];
		optimalsS[0][5] = optimalTable[0][1] + optimalTable[0][2];
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
											.toArray(double[][] :: new);
			int[] numIndex = levelIndex(tOptTable);
			if (numIndex.length == 1 && numIndex[0] != 0) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				if (numIndex[0] == tOptTable.length - 1  && tOptTable[numIndex[0]][2] == maxOrderQuantity) {
					optimalsS[t][0] = tOptTable[numIndex[0]][1] + 1;
					optimalsS[t][1] = tOptTable[numIndex[0]][1] + tOptTable[numIndex[0]][2];
				}
				optimalsS[t][2] = optimalsS[t][0];
				optimalsS[t][3] = optimalsS[t][1];
				optimalsS[t][4] = optimalsS[t][0];
				optimalsS[t][5] = optimalsS[t][1];
			}
			else if (numIndex.length == 1 && numIndex[0] == 0) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0]][1];
				optimalsS[t][2] = tOptTable[numIndex[0]][1];
				optimalsS[t][3] = tOptTable[numIndex[0]][1];
				optimalsS[t][4] = tOptTable[numIndex[0]][1];
				optimalsS[t][5] = tOptTable[numIndex[0]][1];
			}
			else if (numIndex.length == 0) {
				optimalsS[t][0] = tOptTable[tOptTable.length - 1][1];
				optimalsS[t][1] = maxOrderQuantity * 10;
				optimalsS[t][2] = tOptTable[tOptTable.length - 1][1];
				optimalsS[t][3] = maxOrderQuantity * 10;
				optimalsS[t][4] = tOptTable[tOptTable.length - 1][1];
				optimalsS[t][5] = maxOrderQuantity * 10;
			}	
			else if (numIndex.length == 2) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				optimalsS[t][2] = tOptTable[numIndex[1]][1];
				optimalsS[t][3] = tOptTable[numIndex[1] - 1][1] + tOptTable[numIndex[1] - 1][2];
				if (numIndex[1] == tOptTable.length - 1 && tOptTable[numIndex[1]][2] == maxOrderQuantity) {
					optimalsS[t][2] = tOptTable[numIndex[1]][1] + 1;
					optimalsS[t][3] = tOptTable[numIndex[1]][1] + tOptTable[numIndex[1]][2];
				}			
				optimalsS[t][4] = optimalsS[t][2];
				optimalsS[t][5] = optimalsS[t][3];
			}
			else if (numIndex.length == 3) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				optimalsS[t][2] = tOptTable[numIndex[1]][1];
				optimalsS[t][3] = tOptTable[numIndex[1] - 1][1] + tOptTable[numIndex[1] - 1][2];
				optimalsS[t][4] = tOptTable[numIndex[2]][1];
				optimalsS[t][5] = tOptTable[numIndex[2] - 1][1] + tOptTable[numIndex[2] - 1][2];
				if (numIndex[2] == tOptTable.length - 1 && tOptTable[numIndex[2]][2] == maxOrderQuantity) {
					optimalsS[t][4] = tOptTable[numIndex[2]][1] + 1;
					optimalsS[t][5] = tOptTable[numIndex[2]][1] + tOptTable[numIndex[2]][2];
				}
			}
			else {
				int sIndex3 = numIndex[numIndex.length - 1];
				optimalsS[t][4] = tOptTable[sIndex3][1];
				optimalsS[t][5] = tOptTable[sIndex3 - 1][1] + tOptTable[sIndex3 - 1][2];
				int sIndex2 = numIndex[numIndex.length - 2];
				optimalsS[t][2] = tOptTable[sIndex2][1];
				optimalsS[t][3] = tOptTable[sIndex2 - 1][1] + tOptTable[sIndex2 - 1][2];
				int sIndex1 = numIndex[numIndex.length - 3];
				optimalsS[t][0] = tOptTable[sIndex1][1];
				optimalsS[t][1] = minSquare(optimalsS[t][0], sIndex1, tOptTable);
			}
		}
		return optimalsS;
	}
}
