package cash.strongconstraint;

import java.util.ArrayList;
import java.util.Arrays;

import javax.enterprise.inject.Default;

import com.google.common.base.CaseFormat;


/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 14, 2018---1:23:28 PM
 * @description: find s C S for cash flow sdp
 */

public class FindsCS {
	int T;
	double iniCash;

	public FindsCS(int T, double iniCash) {
		this.T = T;
		this.iniCash = iniCash;
	}
	
	public enum FindCCrieria{
		MAX,
		MIN,
		AVG;
	}

	double[][] getsCS(double[][] optimalTable, double minCashRequired, FindCCrieria criteria) {
		double[][] optimalsCS = new double[T][3];
		optimalsCS[0][0] = optimalTable[0][1];
		optimalsCS[0][1] = optimalTable[0][2];
		optimalsCS[0][2] = optimalTable[0][1] + optimalTable[0][3];
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			int recordTimes = 0;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i)
					.map(p -> Arrays.stream(p).toArray()).toArray(double[][]::new);
			optimalsCS[t][2] = minCashRequired;
			optimalsCS[t][1] = 0;
			ArrayList<Double> recordCash = new ArrayList<>();
			double mark_s = 0;
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][3] != 0) {
					if (mark_s == 0) {
						optimalsCS[t][0] = j + 1 < tOptTable.length ? tOptTable[j + 1][1]
																	    : tOptTable[j][1] + 1; // maximum not ordering inventory level as s
						mark_s = 1;
					}
					if (tOptTable[j][1] + tOptTable[j][3] > optimalsCS[t][2]) // a maximum order-up-to level as S
						optimalsCS[t][2] = tOptTable[j][1] + tOptTable[j][3];
					recordTimes = 1;
				}
								
				if (tOptTable[j][3] == 0 && recordTimes == 1) 
					recordCash.add(tOptTable[j][2]);					
			}
			// choose a maximum not ordering cash level as C when ordering quantity is 0
			// or choose an average value
			// or choose a minimum value
			switch (criteria) {
				case MAX:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).max().isPresent() ? recordCash.stream().mapToDouble(p -> p).max().getAsDouble() : minCashRequired;
					break;
				case MIN:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).min().isPresent() ? recordCash.stream().mapToDouble(p -> p).min().getAsDouble() : minCashRequired;
					break;
				default:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).average().isPresent() ? recordCash.stream().mapToDouble(p -> p).average().getAsDouble() : minCashRequired;
					
			}
		}
		System.out.println("(s, C, S) are: " + Arrays.deepToString(optimalsCS));
		return optimalsCS;
	}

	/**
	 * check whether (s, C, S) policy satisfy all states in the optimal table of sdp
	 * 
	 * @param sCS
	 * @param optTable
	 * @param fixOrderCost
	 * @param variCost
	 */
	public int checksBS(double[][] sCS, double[][] optTable, Double minCashRequired, double maxOrderQ, double fixOrderCost, double variCost) {
		int nonOptCount = 0;
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			
			double[][] tOptTable = Arrays.stream(optTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
					.toArray(double[][]::new);
			for (int j = 0; j < tOptTable.length; j++) {
				if (tOptTable[j][1] >= sCS[t][0] && tOptTable[j][3] != 0)
					nonOptCount++;
				if (tOptTable[j][2] <= sCS[t][1] && tOptTable[j][3] != 0)
					nonOptCount++;
				double maxQ = Math.min(sCS[t][2] - tOptTable[j][1], (tOptTable[j][2] - minCashRequired - fixOrderCost) / variCost);
				maxQ = Math.min(maxQ, maxOrderQ);
				if (tOptTable[j][1] < sCS[t][0] && tOptTable[j][2] > sCS[t][1] && tOptTable[j][3] != maxQ)
					nonOptCount++;
			}
		}
		System.out.println("there are " + nonOptCount + " states that not satisfy (s, C, S) ordering property");
		return nonOptCount;
	}
	
}
