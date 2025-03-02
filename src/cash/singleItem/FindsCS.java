package cash.singleItem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.IntStream;

import sdp.inventory.State;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;


/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 14, 2018---1:23:28 PM
 * @description: find s C S for cash flow sdp
 */

public class FindsCS {
	int T;
	double iniCash;
	Distribution[] distributions;
	double fixOrderCost;
	double price;
	double variOrderCost;
	double holdCost;
	double salvageValue;
	
	Map<State, Double> cacheC1Values = new TreeMap<>(); // record C1 values for different initial inventory x
	Map<State, Double> cacheC2Values = new TreeMap<>(); // record C2 values for different initial inventory x
	
	public FindsCS(double iniCash, Distribution[] distributions, double fixOrderCost, double price,
			double variOrderCost, double holdCost, double salvageValue) {
		this.T = distributions.length;
		this.iniCash = iniCash;
		this.distributions = distributions;
		this.fixOrderCost = fixOrderCost;
		this.price = price;
		this.variOrderCost = variOrderCost;
		this.holdCost = holdCost;
		this.salvageValue = salvageValue;
		Comparator<State> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? 0 : -1 : -1;
		this.cacheC1Values = new TreeMap<>(keyComparator);
		this.cacheC2Values = new TreeMap<>(keyComparator);
	}
	
	public enum FindCCrieria{
		MAX,
		MIN,
		AVG,
		XRELATE;
	}
	
	
	// C bound is a big factor causing optimality gaps. Find C through searching in tOptTable, not by computing Ly
	/**
	 * get s, C1(x), C2(x), S
	 */
	double[][] getsC12S(double[][] optimalTable, double overheadCost, FindCCrieria criteria) {
		int M = 10000;
		cacheC1Values.clear();
		cacheC2Values.clear();
		double[][] optimalsCS = new double[T][4];
		optimalsCS[0][0] = optimalTable[0][1] + 1;
		optimalsCS[0][1] = optimalTable[0][2] ;
		optimalsCS[0][2] = M; 
		optimalsCS[0][3] = optimalTable[0][1] + optimalTable[0][3];	
		cacheC1Values.put(new State(1, optimalTable[0][0]), optimalTable[0][1]);
		cacheC2Values.put(new State(1, optimalTable[0][0]), optimalTable[0][2]);
		
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i)
					.map(p -> Arrays.stream(p).toArray()).toArray(double[][]::new);
			
			if (t == T - 1) {
				Distribution distribution = distributions[T - 1];
				optimalsCS[T - 1][3] = distribution.inverseF((price - variOrderCost) / (holdCost  + price - salvageValue));
				optimalsCS[T - 1][2] = M;			
				double S = optimalsCS[T - 1][3];
				for (int j = (int) S; j >= 0; j--) {
					if (Ly(j, t) < Ly(S, t) - fixOrderCost) {
						optimalsCS[T - 1][0] = j + 1;
						break;
					}
				}
				optimalsCS[T - 1][1] = 0; // C default value is 0
				for (int j = (int) S; j >= 0; j--) {
					int jj = 0;
					for (jj = j + 1; jj <= (int) S; jj++) {
						if (Ly(jj,  t) > fixOrderCost + Ly(j, t)) {
							optimalsCS[t][1] = fixOrderCost + variOrderCost * (jj - 1 - j); // C for x = 0 at last period
							cacheC1Values.put(new State(t + 1, j), optimalsCS[t][1]);
							cacheC2Values.put(new State(t + 1, j), optimalsCS[t][2]);
							break;
						}
					}
					if (Ly(S, t) < fixOrderCost) { // choose a large value for C1, since expected profit is too small
						optimalsCS[t][1] = M;
						cacheC1Values.put(new State(t + 1, j), optimalsCS[t][1]);
						cacheC2Values.put(new State(t + 1, j), optimalsCS[t][2]);
					}
				}
				break; // when t = T - 1, no need to compute C and S below
			}				
			
			optimalsCS[t][3] = 0;  // default value for S is 0
			optimalsCS[t][2] = M; // default value for C2 is M
			optimalsCS[t][1] = fixOrderCost + overheadCost; // default value for C1 is K plus overhead cost
			optimalsCS[t][0] = 0;  // default value for s is 0
			ArrayList<Double> recordCash = new ArrayList<>();									
			
			// backward, the first inventory level that starts ordering is s
			boolean sHasRecorded = false; //backward, the first inventory level that starts ordering is s
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][3] != 0) {
					if (sHasRecorded == false) {
						optimalsCS[t][0] = j + 1 < tOptTable.length ? tOptTable[j][1] + 1 
																	    : tOptTable[j][1] + 1; // maximum not ordering inventory level as s
						//optimalsCS[t][0] = tOptTable[j][1] < 1 ? 0 : tOptTable[j][1] + 1; // avoid large s, make s be 0
						
						// choose a most frequent s in this inventory level
//						int zeros = 0;
//						int nonZeros = 1;
//						for (int m = j; m >=0; m--) {
//							if (m < j) {
//								if (tOptTable[m][3] != 0)
//									nonZeros++;
//								else {
//									zeros++;
//								}
//							}
//						}
//						if (zeros>nonZeros)
//							optimalsCS[t][0] = 0;
						sHasRecorded = true;
					}
					if (tOptTable[j][1] + tOptTable[j][3] > optimalsCS[t][3]) //  maximum order-up-to level as S
						optimalsCS[t][3] = tOptTable[j][1] + tOptTable[j][3]; // a template value for S
						//if (tOptTable[j][2] > fixOrderCost + variOrderCost * tOptTable[j][3])
							//recordS.add(tOptTable[j][1] + tOptTable[j][3]); // average order-up-to level as S, is worse than choosing maximum S
					sHasRecorded = true; // 
				}
								
				if (tOptTable[j][3] == 0 && sHasRecorded == true) // for computing C(x)
					recordCash.add(tOptTable[j][2]);					
			}
			
			// choose a maximum not ordering cash level as C when ordering quantity is 0
			// or choose an average value
			// or choose a minimum value
			// or related with x (initial inventory)			
			switch (criteria) {
				case MAX:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).max().isPresent() ? recordCash.stream().mapToDouble(p -> p).max().getAsDouble() : overheadCost;
					break;
				case MIN:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).min().isPresent() ? recordCash.stream().mapToDouble(p -> p).min().getAsDouble() : overheadCost;
					break;
				case AVG:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).average().isPresent() ? recordCash.stream().mapToDouble(p -> p).average().getAsDouble() : overheadCost;
					break;
				case XRELATE:
					// C1
					double markInventory = -0.5;
					boolean CHasRecoded = false;
					for (int j = 0; j < tOptTable.length - 1; j++) {
						if (tOptTable[j][1] < optimalsCS[t][0]) {							
							if (tOptTable[j][3] > 0) {
								if (tOptTable[j][1] > markInventory) {
									optimalsCS[t][1] = tOptTable[j][2] - 1;
									markInventory = tOptTable[j][1];
									CHasRecoded = true;
								}
							}
							else { // if for an initial inventory x, order quantity always zero, then C is a large number
								optimalsCS[t][1] = CHasRecoded == false ? fixOrderCost * 20 : optimalsCS[t][1];
							}
							cacheC1Values.put(new State(t + 1, tOptTable[j][1]), optimalsCS[t][1]);
						}
						else 
							break;
					}
					
					// C2
					markInventory = M;
					CHasRecoded = false;
					for (int j = tOptTable.length - 2; j >= 0; j--) {
						if (tOptTable[j][1] < optimalsCS[t][0]) {							
							if (tOptTable[j][3] > 0 && tOptTable[j + 1][3] == 0 && tOptTable[j][1] == tOptTable[j + 1][1]) 
							{
								if (tOptTable[j][1] < markInventory ) {
									optimalsCS[t][2] = tOptTable[j][2] + 1;
									markInventory = tOptTable[j][1];
									CHasRecoded = true;
									cacheC2Values.put(new State(t + 1, tOptTable[j][1]), optimalsCS[t][2]);
								}
								//cacheC2Values.put(new State(t + 1, tOptTable[j][1]), optimalsCS[t][2]);
							}
							
						}
					}
					break;
			}
					
			// find a most frequent S, sometimes when cash is not enough, S bound can also affect gaps much
			//Comparator<Map.Entry<Double, Integer>> comparator = (o1, o2) -> o1.getValue() > o2.getValue() ? 1
			//																: o1.getValue() == o2.getValue() ?
			//																  o1.getKey() > o2.getKey() ? 1 : -1 : -1;
			Map<Double, Integer> recordS = new HashMap<>();
			double S = 0;
			boolean lastOrderFullCapacity = false;
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][1] < optimalsCS[t][0] && tOptTable[j][3] > 0) {
					int maxQ = (int) Math.max(0, (tOptTable[j][2] - overheadCost - fixOrderCost) / variOrderCost);
					if (tOptTable[j][2] >= fixOrderCost + variOrderCost * tOptTable[j][3] + overheadCost) {														
							if (recordS.size() != 0) {
								S = (double) recordS.keySet().toArray()[recordS.size() - 1];
								if (tOptTable[j][1] + tOptTable[j][3] > S + 0.1) { // order at full capacity, real S may not be increasing
									if (lastOrderFullCapacity == true) {
										int num = recordS.get(S) + 1;          
										recordS.remove(S);
										S = tOptTable[j][1] + tOptTable[j][3];									
										recordS.putIfAbsent(S, num);
									}
									else {
//										S = tOptTable[j][1] + tOptTable[j][3];									
//										recordS.putIfAbsent(S, 1);
//										lastOrderFullCapacity = tOptTable[j][3] > maxQ - 0.1 ? true : false;
										int num = recordS.get(S) + 1; 
										recordS.remove(S);
										S = tOptTable[j][1] + tOptTable[j][3];	// change 1  to num
										recordS.putIfAbsent(S, num);
										lastOrderFullCapacity = tOptTable[j][3] > maxQ - 0.1 ? true : false;		
									}									
								}
								else {
									if (tOptTable[j][3] < maxQ - 0.1 && tOptTable[j][1] + tOptTable[j][3] < S - 0.1) { //new S
										
										recordS.putIfAbsent(tOptTable[j][1] + tOptTable[j][3], 1);
										lastOrderFullCapacity = false;
									} 										
									else // order at full capacity or same S
										recordS.replace(S, recordS.get(S) + 1);
								}
							}
							else {
								S = tOptTable[j][1] + tOptTable[j][3];
								recordS.putIfAbsent(S, 1);
								if (tOptTable[j][3] > maxQ - 0.1)
									lastOrderFullCapacity = true;
							}
					}	
				}
			}			
			if (recordS.size() != 0) 
				optimalsCS[t][3] = recordS.entrySet().stream()
				.max((o1,o2)-> o1.getValue() > o2.getValue() ? 1 :
					o1.getValue() == o2.getValue() ?  
							o1.getKey() > o2.getKey() ? 1 : -1 : -1).get().getKey();
			if (recordS.size() == 0 && optimalsCS[t][0] != 0)
				optimalsCS[t][3] = M; // a large number when ordering quantity is always full capacity
															
					//double meands = t == T - 1 ? meanD[t] : meanD[t] + meanD[t + 1]; // a heuristic step
//					double meands = meanD[t];
//					Distribution distribution = new PoissonDist(meands);
//					double S = distribution.inverseF((price - variOrderCost) / (holdCost  + price));
//					for (int j = 0; j <= (int) S; j++) {
//						int jj = 0;
//						for (jj = j + 1; jj <= (int) S; jj++) {
//							if (Ly(jj,  t) > fixOrderCost + Ly(j, t)) {
//								optimalsCS[t][1] = fixOrderCost + variOrderCost * (jj - 1 - j); // C for x = 0 at last
//								cacheC1Values.put(new State(t + 1, j), optimalsCS[t][1]);
//								break;
//							}				
//							
//						}
//						if (Ly(S, t) < fixOrderCost) { // choose a large value for C, since expected profit is too small
//							optimalsCS[t][1] = fixOrderCost * 20;
//							cacheC1Values.put(new State(t + 1, j), optimalsCS[t][1]);
//						}
//					}									
		}
		
  		System.out.println("(s, C1, C2, S) are: " + Arrays.deepToString(optimalsCS));
		return optimalsCS;
	}
	
	
	/**
	 * get s, C1(x), S, 
	 * can choose most frequent S, or max S when without inventory holding cost
	 */
	
	double[][] getsCS(double[][] optimalTable, double overheadCost, FindCCrieria criteria) {
		int M = 10000;
		cacheC1Values.clear();
		cacheC2Values.clear();
		double[][] optimalsCS = new double[T][4];
		optimalsCS[0][0] = optimalTable[0][1] + 1;
		optimalsCS[0][1] = 0; // optimalTable[0][2] ;
		optimalsCS[0][2] = M; 
		optimalsCS[0][3] = optimalTable[0][1] + optimalTable[0][3];	
		cacheC1Values.put(new State(1, optimalTable[0][0]), optimalTable[0][1]);
		cacheC2Values.put(new State(1, optimalTable[0][0]), optimalTable[0][2]);
		
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i)
					.map(p -> Arrays.stream(p).toArray()).toArray(double[][]::new);
			
			if (t == T - 1) {
				Distribution distribution = distributions[T - 1];
				optimalsCS[T - 1][3] = variOrderCost != salvageValue ? distribution.inverseF((price - variOrderCost) / (holdCost  + price - salvageValue))
																	 : distribution.inverseF(0.999);
				optimalsCS[T - 1][2] = M;			
				double S = optimalsCS[T - 1][3];
				optimalsCS[T - 1][0] = 0; // s default value is 0
				for (int j = (int) S; j >= 0; j--) {
					if (Ly(j, t) < Ly(S, t) - fixOrderCost) {
						optimalsCS[T - 1][0] = j + 1;
						break;
					}
				}
				optimalsCS[T - 1][1] = 0; // C default value is M
				for (int j = (int) S; j >= 0; j--) {
					int jj = 0;
					for (jj = j + 1; jj <= (int) S; jj++) {
						if (Ly(jj,  t) > fixOrderCost + Ly(j, t)) {
							optimalsCS[t][1] = fixOrderCost + variOrderCost * (jj - 1 - j); // C for x = 0 at last period
							cacheC1Values.put(new State(t + 1, j), optimalsCS[t][1]);
							break;
						}
					}
					if (Ly(S, t) < fixOrderCost) { // choose a large value for C1, since expected profit is too small
						optimalsCS[t][1] = M;
						cacheC1Values.put(new State(t + 1, j), optimalsCS[t][1]);
					}
				}
				break; // when t = T - 1, no need to compute C and S below
			}				
			
			optimalsCS[t][3] = 0;  // default value for S is 0
			optimalsCS[t][2] = M; // default value for C2 is M
			optimalsCS[t][1] = fixOrderCost + overheadCost; // default value for C1 is K
			optimalsCS[t][0] = 0;  // default value for s is 0
			ArrayList<Double> recordCash = new ArrayList<>();									
			
			// backward, the first inventory level that starts ordering is s - 1
			boolean sHasRecorded = false; //backward, the first inventory level that starts ordering is s + 1
			                              // or s is the minimum initial inventory level that results in not ordering
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][3] != 0) {
					if (sHasRecorded == false) {
						optimalsCS[t][0] = j + 1 < tOptTable.length ? tOptTable[j][1] + 1 
																	    : tOptTable[j][1] + 1; // maximum not ordering inventory level as s
						
						// choose a most frequent s in this inventory level
//						int zeros = 0;
//						int nonZeros = 1;
//						for (int m = j; m >=0; m--) {
//							if (m < j) {
//								if (tOptTable[m][3] != 0)
//									nonZeros++;
//								else {
//									zeros++;
//								}
//							}
//						}
//						if (zeros>nonZeros)
//							optimalsCS[t][0] = 0;
						
						sHasRecorded = true;
					}
//					optimalsCS[2][0] =0;
//					if (tOptTable[j][1] + tOptTable[j][3] > optimalsCS[t][3]) //  maximum order-up-to level as S
//						optimalsCS[t][3] = tOptTable[j][1] + tOptTable[j][3];
					// the following codes to compute S when h is 0. It can also be a temperate value of S
					double demandSum = IntStream.range(t, T).mapToObj(k -> distributions[k].getMean())
							.reduce(0.0, (x, y) -> x.doubleValue() + y.doubleValue());
					Distribution distribution = new PoissonDist(demandSum);
					optimalsCS[t][3]  = distribution.inverseF((price - variOrderCost) / (price - salvageValue));
						//if (tOptTable[j][2] > fixOrderCost + variOrderCost * tOptTable[j][3])
							//recordS.add(tOptTable[j][1] + tOptTable[j][3]); // average order-up-to level as S, is worse than choosing maximum S
					sHasRecorded = true; 
				}
								
				if (tOptTable[j][3] == 0 && sHasRecorded == true) 
					recordCash.add(tOptTable[j][2]);					
			}
			
			// choose a maximum not ordering cash level as C when ordering quantity is 0
			// or choose an average value
			// or choose a minimum value
			// or related with x (initial inventory)			
			switch (criteria) {
				case MAX:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).max().isPresent() ? recordCash.stream().mapToDouble(p -> p).max().getAsDouble() : overheadCost;
					break;
				case MIN:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).min().isPresent() ? recordCash.stream().mapToDouble(p -> p).min().getAsDouble() : overheadCost;
					break;
				case AVG:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).average().isPresent() ? recordCash.stream().mapToDouble(p -> p).average().getAsDouble() : overheadCost;
					break;
				case XRELATE:
					// C1
					double markInventory = -0.5;
					boolean CHasRecoded = false;
					for (int j = 0; j < tOptTable.length; j++) { // the length of j
						if (tOptTable[j][1] < optimalsCS[t][0]) {							
							if (tOptTable[j][3] > 0) {
								if (tOptTable[j][1] > markInventory) {
									optimalsCS[t][1] = tOptTable[j][2] - 1;
									markInventory = tOptTable[j][1];
									CHasRecoded = true;
								}
							}
							else { // if for an initial inventory x, order quantity always zero, then C is a large number
								optimalsCS[t][1] = CHasRecoded == false ? fixOrderCost * 20 : optimalsCS[t][1];
							}
							cacheC1Values.put(new State(t + 1, tOptTable[j][1]), optimalsCS[t][1]);
						}
						else 
							break;
					}
					break;
			}
					
			// find a most frequent S, sometimes when cash is not enough, S bound can also affect gaps much
			//Comparator<Map.Entry<Double, Integer>> comparator = (o1, o2) -> o1.getValue() > o2.getValue() ? 1
			//																: o1.getValue() == o2.getValue() ?
			//																  o1.getKey() > o2.getKey() ? 1 : -1 : -1;
			// this is used when having positive holding cost h
			Map<Double, Integer> recordS = new HashMap<>();
			double S = 0;
			boolean lastOrderFullCapacity = false;
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][1] < optimalsCS[t][0] && tOptTable[j][3] > 0) {
					int maxQ = (int) Math.max(0, (tOptTable[j][2] - overheadCost - fixOrderCost) / variOrderCost);
					if (tOptTable[j][2] >= fixOrderCost + variOrderCost * tOptTable[j][3] + overheadCost) {														
						if (recordS.size() != 0) {
							S = (double) recordS.keySet().toArray()[recordS.size() - 1];
							if (tOptTable[j][1] + tOptTable[j][3] >= S) { // order at full capacity, real S may not be increasing
								if (lastOrderFullCapacity == true) {
									int num = recordS.get(S) + 1;          
									recordS.remove(S);
									S = tOptTable[j][1] + tOptTable[j][3];									
									recordS.putIfAbsent(S, num);
								}
								else {
									int num = recordS.get(S) + 1; 
									recordS.remove(S);
									S = tOptTable[j][1] + tOptTable[j][3];	// change 1  to num
									recordS.putIfAbsent(S, num);
									lastOrderFullCapacity = tOptTable[j][3] > maxQ - 0.1 ? true : false;										
								}									
							}
							else {
								if (tOptTable[j][3] < maxQ - 0.1 && tOptTable[j][1] + tOptTable[j][3] < S - 0.1) { //new S

									recordS.putIfAbsent(tOptTable[j][1] + tOptTable[j][3], 1);
									lastOrderFullCapacity = false;
								} 										
								else // order at full capacity or same S
									recordS.replace(S, recordS.get(S) + 1);
							}
						}
						else {
							S = tOptTable[j][1] + tOptTable[j][3];
							recordS.putIfAbsent(S, 1);
							if (tOptTable[j][3] > maxQ - 0.1)
								lastOrderFullCapacity = true;
						}
					}	
				}
			}			
			if (recordS.size() != 0) 
				optimalsCS[t][3] = recordS.entrySet().stream()
				.max((o1,o2)-> o1.getValue() > o2.getValue() ? 1 :
					o1.getValue() == o2.getValue() ?  
							o1.getKey() > o2.getKey() ? 1 : -1 : -1).get().getKey();
			if (recordS.size() == 0 && optimalsCS[t][0] != 0)
				optimalsCS[t][3] = M; // a large number when ordering quantity is always full capacity
		}		
		
		
		double[][] sCS = new double[T][3];
		for (int t = 0; t < T; t++) {
			sCS[t][0] = optimalsCS[t][0];
			sCS[t][1] = optimalsCS[t][1];
			sCS[t][2] = optimalsCS[t][3];
		}	
  		System.out.println("(s, C, S) are: " + Arrays.deepToString(sCS));
		return sCS;
	}

	
	/**
	 * @param optimalTable
	 * @param overheadCost
	 * @param criteria
	 * @return get the number of S levels and the frequency of the most frequent one in each period 
	 * @date: May 2, 2020, 10:55:39 AM 
	 */
	double[][] getMaxSFrequency(double[][] optimalTable, double overheadCost, FindCCrieria criteria){
		double[][] numFrequency = new double[T][2];
		numFrequency[0][0] = 1; numFrequency[0][1] = 100;
		numFrequency[T - 1][0] = 1; numFrequency[T - 1][1] = 100;
		
		double[][] optimalsCS = new double[T][4];
		optimalsCS[0][0] = optimalTable[0][1] + 1;
		optimalsCS[0][1] = 0; // optimalTable[0][2] ;
		int M = 10000;
		optimalsCS[0][2] = M; 
		optimalsCS[0][3] = optimalTable[0][1] + optimalTable[0][3];	
		cacheC1Values.put(new State(1, optimalTable[0][0]), optimalTable[0][1]);
		cacheC2Values.put(new State(1, optimalTable[0][0]), optimalTable[0][2]);
		
		for (int t = 1; t < T - 1; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i)
					.map(p -> Arrays.stream(p).toArray()).toArray(double[][]::new);
			
			
			optimalsCS[t][3] = 0;  // default value for S is 0
			optimalsCS[t][2] = M; // default value for C2 is M
			optimalsCS[t][1] = fixOrderCost + overheadCost; // default value for C1 is K
			optimalsCS[t][0] = 0;  // default value for s is 0
			ArrayList<Double> recordCash = new ArrayList<>();									
			
			// backward, the first inventory level that starts ordering is s - 1
			boolean sHasRecorded = false; //backward, the first inventory level that starts ordering is s
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][3] != 0) {
					if (sHasRecorded == false) {
						optimalsCS[t][0] = j + 1 < tOptTable.length ? tOptTable[j][1] + 1 
																	    : tOptTable[j][1] + 1; // maximum not ordering inventory level as s
						sHasRecorded = true;
					}
					if (tOptTable[j][1] + tOptTable[j][3] > optimalsCS[t][3]) //  maximum order-up-to level as S
						optimalsCS[t][3] = tOptTable[j][1] + tOptTable[j][3]; // this S is template value
						//if (tOptTable[j][2] > fixOrderCost + variOrderCost * tOptTable[j][3])
							//recordS.add(tOptTable[j][1] + tOptTable[j][3]); // average order-up-to level as S, is worse than choosing maximum S
					sHasRecorded = true; // 
				}
								
				if (tOptTable[j][3] == 0 && sHasRecorded == true) 
					recordCash.add(tOptTable[j][2]);					
			}
			
			// choose a maximum not ordering cash level as C when ordering quantity is 0
			// or choose an average value
			// or choose a minimum value
			// or related with x (initial inventory)			
			switch (criteria) {
				case MAX:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).max().isPresent() ? recordCash.stream().mapToDouble(p -> p).max().getAsDouble() : overheadCost;
					break;
				case MIN:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).min().isPresent() ? recordCash.stream().mapToDouble(p -> p).min().getAsDouble() : overheadCost;
					break;
				case AVG:
					optimalsCS[t][1] = recordCash.stream().mapToDouble(p -> p).average().isPresent() ? recordCash.stream().mapToDouble(p -> p).average().getAsDouble() : overheadCost;
					break;
				case XRELATE:
					// C1
					double markInventory = -0.5;
					boolean CHasRecoded = false;
					for (int j = 0; j < tOptTable.length; j++) { // the length of j
						if (tOptTable[j][1] < optimalsCS[t][0]) {							
							if (tOptTable[j][3] > 0) {
								if (tOptTable[j][1] > markInventory) {
									optimalsCS[t][1] = tOptTable[j][2] - 1;
									markInventory = tOptTable[j][1];
									CHasRecoded = true;
								}
							}
							else { // if for an initial inventory x, order quantity always zero, then C is a large number
								optimalsCS[t][1] = CHasRecoded == false ? fixOrderCost * 20 : optimalsCS[t][1];
							}
							cacheC1Values.put(new State(t + 1, tOptTable[j][1]), optimalsCS[t][1]);
						}
						else 
							break;
					}
					break;
			}
			
			Map<Double, Integer> recordS = new HashMap<>();
			double S = 0;
			boolean lastOrderFullCapacity = false;
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][1] < optimalsCS[t][0] && tOptTable[j][3] > 0) {
					int maxQ = (int) Math.max(0, (tOptTable[j][2] - overheadCost - fixOrderCost) / variOrderCost);
					if (tOptTable[j][2] >= fixOrderCost + variOrderCost * tOptTable[j][3] + overheadCost) {														
						if (recordS.size() != 0) {
							S = (double) recordS.keySet().toArray()[recordS.size() - 1];
							if (tOptTable[j][1] + tOptTable[j][3] >= S) { // order at full capacity, real S may not be increasing
								if (lastOrderFullCapacity == true) {
									int num = recordS.get(S) + 1;          
									recordS.remove(S);
									S = tOptTable[j][1] + tOptTable[j][3];									
									recordS.putIfAbsent(S, num);
								}
								else {
									int num = recordS.get(S) + 1; 
									recordS.remove(S);
									S = tOptTable[j][1] + tOptTable[j][3];	// change 1  to num
									recordS.putIfAbsent(S, num);
									lastOrderFullCapacity = tOptTable[j][3] > maxQ - 0.1 ? true : false;										
								}									
							}
							else {
								// new S
								if (tOptTable[j][3] > 0.1 && tOptTable[j][3] < maxQ - 0.1 && tOptTable[j][1] + tOptTable[j][3] < S - 0.1) { //new S

									recordS.putIfAbsent(tOptTable[j][1] + tOptTable[j][3], 1);
									lastOrderFullCapacity = false;
								} 										
								else // order at full capacity or same S, deem as same S
									recordS.replace(S, recordS.get(S) + 1); // S and its frequency
							}
						}
						else {
							S = tOptTable[j][1] + tOptTable[j][3];
							recordS.putIfAbsent(S, 1);
							if (tOptTable[j][3] > maxQ - 0.1)
								lastOrderFullCapacity = true;
						}
					}	
				}
			}			
			if (recordS.size() != 0) {
				Iterator<Map.Entry<Double, Integer>> iterator = recordS.entrySet().iterator();
				numFrequency[t][0] = recordS.size();
				double maxS = 0; double maxFre = 0; double totalNum = 0;
				while (iterator.hasNext()) {
					Map.Entry<Double, Integer> entry = iterator.next();
					if (entry.getKey() > maxS) {
						maxS = entry.getKey();
						maxFre =  entry.getValue();
					}
					totalNum += entry.getValue();
				}
				numFrequency[t][1] = maxFre * 100 / totalNum;
				if (recordS.size() > 1) // without holding cost, there will only be one S
					System.out.println(recordS.toString());
			}
			if (recordS.size() == 0 && optimalsCS[t][0] != 0) { // when ordering quantity is always full capacity
				numFrequency[t][0] = 1; numFrequency[t][1] = 100;
			}		
		}	
		
		return numFrequency;
	}
	
		
	/**
	 * compute L(y) in a single period
	 * @param y : order-up-to level y,
	 * @param t : period t - 1
	 */
	double Ly(double y, int t) {
		Distribution distribution = distributions[t];
		
		double meanI = 0;
		for (int i = 0; i < y; i++)
			meanI += (y - i) * (distribution.cdf(i + 0.5) - distribution.cdf(i - 0.5));
		
		double Ly = 0;
		if (t == T - 1)
			Ly = (price - variOrderCost) * y- (price + holdCost - salvageValue) * meanI;
		else
			Ly = (price - variOrderCost) * y- (price + holdCost) * meanI;

		return Ly;
	}
	
	/**
	 * 
	 * @return optimal C values of cacheC1Values
	 */
	public double[][] getOptCTable(){
		Iterator<Map.Entry<State, Double>> iterator = cacheC1Values.entrySet().iterator();
		double[][] arr = new double[cacheC1Values.size()][2];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<State, Double> entry = iterator.next();
			arr[i++] =new double[]{entry.getKey().getPeriod(), entry.getKey().getIniInventory(), entry.getValue()};
		}
		return arr;
	}
	
	/**
	 * check whether (s, C, S) policy satisfy all states in the optimal table of sdp;
	 * when no inventory holding cost, (s, C, S) may not be optimal, because sometimes may be exist C2,
	 * or S is not always follows, S may be x during the states
	 * 
	 * @param sCS
	 * @param optTable
	 * @param fixOrderCost
	 * @param variCost
	 */
	public int checksCS(double[][] sCS, double[][] optTable, Double overheadCost, double maxOrderQ, double fixOrderCost, double variCost) {
		int nonOptCount = 0;
		for (int t = 1; t < T; t++) {
			final int i = t + 1;			
			double[][] tOptTable = Arrays.stream(optTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
					.toArray(double[][]::new);
			for (int j = 0; j < tOptTable.length; j++) {
				if (tOptTable[j][1] < sCS[t][0]) {
					if (cacheC1Values.get(new State(i, tOptTable[j][1])) == null)
						sCS[t][1] = 0;
					else
						sCS[t][1] = cacheC1Values.get(new State(i, tOptTable[j][1]));
				}			
				if (tOptTable[j][1] >= sCS[t][0] && tOptTable[j][3] != 0) // check s
					nonOptCount++;
				if (tOptTable[j][2] <= sCS[t][1] && tOptTable[j][3] != 0) // check C
					nonOptCount++;
				double maxQ = (int) Math.min(sCS[t][2] - tOptTable[j][1], (tOptTable[j][2] - overheadCost - fixOrderCost) / variCost);
				maxQ = Math.min(maxQ, maxOrderQ);	
				if (tOptTable[j][1] < sCS[t][0] && tOptTable[j][2] > sCS[t][1] && tOptTable[j][3] != maxQ)
					nonOptCount++;
			}
		}
		long totalStateNum = optTable.length;
		System.out.println("there are " + nonOptCount + "/" + totalStateNum + " states that not satisfy (s, C, S) ordering property");
		return nonOptCount;
	}
	
	/**
	 * check whether (s, C1, C2, S) policy satisfy all states in the optimal table of sdp
	 * 
	 * @param sC1C2S
	 * @param optTable
	 * @param fixOrderCost
	 * @param variCost
	 */
	public int checksC12S(double[][] sCS, double[][] optTable, Double overheadCost, double maxOrderQ, double fixOrderCost, double variCost) {
		int nonOptCount = 0;
		for (int t = 1; t < T; t++) {
			final int i = t + 1;			
			double[][] tOptTable = Arrays.stream(optTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
					.toArray(double[][]::new);
			for (int j = 0; j < tOptTable.length; j++) {
				if (tOptTable[j][1] < sCS[t][0]) {
					if (cacheC1Values.get(new State(i, tOptTable[j][1])) == null)
						sCS[t][1] = 0;
					else
						sCS[t][1] = cacheC1Values.get(new State(i, tOptTable[j][1]));
					if (cacheC2Values.get(new State(i, tOptTable[j][1])) == null)
						sCS[t][2] = 0;
					else
						sCS[t][2] = cacheC2Values.get(new State(i, tOptTable[j][1]));
				}			
				if (tOptTable[j][1] >= sCS[t][0] && tOptTable[j][3] != 0)
					nonOptCount++;
				if (tOptTable[j][2] <= sCS[t][1] && tOptTable[j][3] != 0)
					nonOptCount++;
				if (tOptTable[j][2] >= sCS[t][2] && tOptTable[j][3] != 0)
					nonOptCount++;
				double maxQ = (int) Math.min(sCS[t][3] - tOptTable[j][1], (tOptTable[j][2] - overheadCost - fixOrderCost) / variCost);
				maxQ = Math.min(maxQ, maxOrderQ);	
				
				if (tOptTable[j][1] < sCS[t][0] && tOptTable[j][2] > sCS[t][1] && tOptTable[j][2] < sCS[t][2] && tOptTable[j][3] != maxQ)
					nonOptCount++;
			}
		}
		long totalStateNum = optTable.length;
		System.out.println("there are " + nonOptCount + "/" + totalStateNum + " states that not satisfy (s, C1, C2, S) ordering property");
		return nonOptCount;
	}
	

	
}
