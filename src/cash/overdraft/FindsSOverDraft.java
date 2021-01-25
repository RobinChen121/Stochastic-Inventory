package cash.overdraft;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 *@author: Zhen Chen
 *@email: 15011074486@163.com
 *@date: Jul 15, 2018---10:36:21 AM
 *@description:  find approximate s, S for overdraft sdp problem
 */

public class FindsSOverDraft {
	int T;
	double iniCash;

	public FindsSOverDraft(int T, double iniCash) {
		this.T = T;
		this.iniCash = iniCash;
	}

	
	/**
	 * find s, S values for overdraft problem
	 * @param optimalTable
	 * @return
	 */
	public double[][] getsS(double[][] optimalTable) {
		double[][] optimalsS = new double[T][2];
		optimalsS[0][0] = optimalTable[0][1];
		optimalsS[0][1] = optimalTable[0][1] + optimalTable[0][3];
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] ==
					i).map(p -> Arrays.stream(p).toArray())
					.toArray(double[][] :: new);

			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][3] != 0) {
					optimalsS[t][0] = tOptTable[j+1][1];
					optimalsS[t][1] = tOptTable[j][1] + tOptTable[j][3];
					break;
				}
			}
		}
		System.out.println("(s, S) are: " + Arrays.deepToString(optimalsS));
		return optimalsS;
	}
	
	/**
	 * find s, C, S for overdraft with limit
	 * @param optimalTable
	 * @return
	 */
	public double[][] getsCS(double[][] optimalTable) {
		double[][] optimalsCS = new double[T][3];
		double minCashUsed = -500;
		ArrayList<Double> recordC = new ArrayList<>(); 
		optimalsCS[0][0] = optimalTable[0][1];
		optimalsCS[0][1] = optimalTable[0][2];
		optimalsCS[0][2] = optimalTable[0][1] + optimalTable[0][3];
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			int recordsTime = 0;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] ==
					i).map(p -> Arrays.stream(p).toArray())
					.toArray(double[][] :: new);

			optimalsCS[t][2] = 0;
			optimalsCS[t][1] = minCashUsed; 
			double mark_s = 0;
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][3] != 0) {
					if (mark_s == 0) {
						optimalsCS[t][0] = tOptTable[j + 1][1]; // maximum not ordering inventory level as s
						mark_s = 1;
					}
					if (tOptTable[j][2] + tOptTable[j][3] > optimalsCS[t][2])  
						optimalsCS[t][2] = tOptTable[j][1] + tOptTable[j][3]; // a maximum order-up-to level as S
					recordsTime = 1;
				}
				
				// choose a maximum not ordering cash level as C when ordering quantity is 0
				// sometimes choose an average value
				if (tOptTable[j][3] == 0 && recordsTime == 1) { 
					recordC.add(tOptTable[j][2]);					
					if (tOptTable[j][2] > optimalsCS[t][1]) // choose maximum value
						optimalsCS[t][1] = tOptTable[j][2];
				}
			}
			//optimalsCS[t][1] = recordC.stream().mapToDouble(s -> s.doubleValue()).sum()/recordC.size();
		}
		System.out.println("(s, C, S) are: " + Arrays.deepToString(optimalsCS));
		return optimalsCS;
	}
	
	
	
	/**
	 * find s, C, S1, S2 for overdraft without limit,
	 * when x < s, order to S1 when cash <= C, order to S2 when cash > C.
	 * 
	 * @param optimalTable
	 * @return s, C, S1, S2
	 */
	public double[][] getsCS1S2(double[][] optimalTable) {
		double[][] optimalsCS1S2 = new double[T][4];
		double minCashUsed = -500;
		optimalsCS1S2[0][0] = optimalTable[0][1]; // s
		optimalsCS1S2[0][1] = optimalTable[0][2]; // C
		optimalsCS1S2[0][2] = optimalTable[0][1] + optimalTable[0][3]; // S1
		optimalsCS1S2[0][3] = optimalTable[0][1] + optimalTable[0][3]; // S2
		

		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] ==
					i).map(p -> Arrays.stream(p).toArray()).toArray(double[][] :: new);

			optimalsCS1S2[t][2] = 0; // S1
			optimalsCS1S2[t][3] = 0; // S2
			optimalsCS1S2[t][1] = minCashUsed; 
			int sIndex = 0;
			
			// compute S frequency for C and S1, S2
			Map<Integer, Integer> recordS = new HashMap<>(); // TreeMap always sorted by keys
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][3] != 0) {
					optimalsCS1S2[t][0] = tOptTable[j + 1][1]; // maximum not ordering inventory level as s
					sIndex = j;
					for (int m = sIndex; m >= 0; m--) {
						if (tOptTable[m][3] != 0) { // whether taking zero Q into account
							int S = (int) (tOptTable[m][1] + tOptTable[m][3]);	
							if (recordS.containsKey(S))
								recordS.replace(S, recordS.get(S) + 1);
							else 
								recordS.putIfAbsent(S, 1);
						}
						
					}	
				break;
				}
			}
			
			// sorting for HashMap by value in descending order
			Comparator<Map.Entry<Integer, Integer>> comparator = (o1, o2) -> o1.getValue() < o2.getValue() ? 1
																			: o1.getValue() == o2.getValue() ?
																			  o1.getKey() < o2.getKey() ? 1 : -1 : -1;
			Map<Integer, Integer> sortedRecordS = recordS.entrySet().stream().sorted(comparator)
					.collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                            (e1, e2) -> e1, LinkedHashMap::new));;
																		  
			// get first two keys
			double S2 = (Integer)sortedRecordS.keySet().toArray()[0];
			double S1 = (Integer)sortedRecordS.keySet().toArray()[1];
			optimalsCS1S2[t][2] = S1; // S1 < S2
			optimalsCS1S2[t][3] = S2; // S2
			
			// find C
			for (int m = sIndex; m >= 0; m--) {
				if (tOptTable[m][1] + tOptTable[m][3] <= S1 + 0.1 && tOptTable[m][1] + tOptTable[m][3] >= S1 - 0.1) {
					optimalsCS1S2[t][1] = tOptTable[m][2];
					break;
				} 
					
			}
		}
		System.out.println("(s, C, S) are: " + Arrays.deepToString(optimalsCS1S2));
		return optimalsCS1S2;
	}
	

}
