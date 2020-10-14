package sdp.cash;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 13, 2018---4:40:05 PM
*@description:  recursion class for cash flow problem
*/

public class CashRecursion {
	Map<CashState, Double> cacheActions = new TreeMap<>();	
	Map<CashState, Double> cacheValues = new TreeMap<>();
	
	double[][][] pmf;	
	OptDirection optDirection;	
	Function<CashState, double[]> getFeasibleActions;
	StateTransitionFunction<CashState, Double, Double, CashState> stateTransition;
	ImmediateValueFunction<CashState, Double, Double, Double> immediateValue;
	double discountFactor;
	
	public enum OptDirection{
		MIN,
		MAX
	}
	
	public CashRecursion(OptDirection optDirection, double[][][] pmf, 
			         Function<CashState, double[]> getFeasibleAction,
			         StateTransitionFunction<CashState, Double, Double, CashState> stateTransition,
			         ImmediateValueFunction<CashState, Double, Double, Double> immediateValue, 
			         double discountFactor) {
		this.optDirection = optDirection;
		this.pmf = pmf;
		this.getFeasibleActions = getFeasibleAction;
		this.stateTransition = stateTransition ;
		this.immediateValue = immediateValue;
		Comparator<CashState> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? o1.iniCash > o2.iniCash  ? 1 :
					o1.iniCash == o2.iniCash ? 0 : -1 : -1 : -1;
		this.cacheActions = new TreeMap<>(keyComparator);
		this.cacheValues = new TreeMap<>(keyComparator);
		this.discountFactor = discountFactor;
	}
		
	public StateTransitionFunction<CashState, Double, Double, CashState> getStateTransitionFunction(){
		return stateTransition;
	}
	
	public ImmediateValueFunction<CashState, Double, Double, Double> getImmediateValueFunction(){
		return immediateValue;
	} 
	
	/**
	 * set a tree map for finding s B S 
	 */
	public void setTreeMapCacheAction() {
		// cacheActions is a sorted map
		Comparator<CashState> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? o1.iniCash > o2.iniCash  ? 1 :
					o1.iniCash == o2.iniCash ? 0 : -1 : -1 : -1;
		cacheActions = new TreeMap<>(keyComparator);
	}
	
	
	public double getExpectedValue(CashState initialState) {
		return this.cacheValues.computeIfAbsent(initialState, s -> {			
//			double val = Arrays.stream(getFeasibleActions.apply(s))
//					//.parallel() // whether using parallel computation, there is error now
//					.map(orderQty -> Arrays.stream(pmf[s.getPeriod() - 1])
//					.mapToDouble(p -> p[1] * immediateValue.apply(s, orderQty, p[0])
//							+ (s.getPeriod() < pmf.length ? p[1] * getExpectedValue(stateTransition.apply(s, orderQty, p[0])) : 0))
//					.sum())
//					.reduce((x, y) -> optDirection == OptDirection.MIN ? x > y ? y : x   // represent min max
//							                                          : x > y ? x :y)
//					.getAsDouble();
//			
//			double bestOrderQty = Arrays.stream(getFeasibleActions.apply(s)).filter(orderQty -> Arrays
//					.stream(pmf[s.getPeriod() - 1])
//					//.parallel() // whether using parallel computation
//					.mapToDouble(p -> p[1] * immediateValue.apply(s, orderQty, p[0])
//							+ (s.getPeriod() < pmf.length ? p[1] * getExpectedValue(stateTransition.apply(s, orderQty, p[0])) : 0))
//					.sum() == val).findAny().getAsDouble();			
			
			double[] feasibleActions = getFeasibleActions.apply(initialState);
			if (initialState.getPeriod() > 1)
				feasibleActions = new double[] {0};
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] QValues = new double[feasibleActions.length];
			double val = optDirection == OptDirection.MIN ? Double.MAX_VALUE
														  : -Double.MAX_VALUE;
		
			double bestOrderQty = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderQty = feasibleActions[i];
				
//				if (s.getPeriod() == 1) { // for debugging
//					orderQty = 0;
//				}
				
				double thisQValue = 0;								
				for (int j = 0; j < dAndP.length; j++) {
					double thisValue = immediateValue.apply(s, orderQty, dAndP[j][0]);
					thisQValue += dAndP[j][1] * immediateValue.apply(s, orderQty, dAndP[j][0]);
					if (s.getPeriod() < pmf.length) {
						CashState newState = stateTransition.apply(s, orderQty, dAndP[j][0]);
						thisQValue += dAndP[j][1] * discountFactor * getExpectedValue(newState);
					}
				}
				QValues[i] = thisQValue;
				if (optDirection == OptDirection.MIN) {
					if (QValues[i] < val) {
						val = QValues[i];
						bestOrderQty = orderQty;
					}
				}
				else {
					if (QValues[i] > val) {
						val = QValues[i];
						bestOrderQty = orderQty;
					}
				}
			}
			try {
			this.cacheActions.putIfAbsent(s, bestOrderQty);
			}
			catch (Exception e) {
				System.out.println("error");
			}
			return val;
		});
	}
	
	public double getAction(CashState state) {
		return cacheActions.get(state);
	}
	
	public Map<CashState, Double> getCacheActions() {
		return cacheActions;
	}
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(){
		Iterator<Map.Entry<CashState, Double>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][3];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<CashState, Double> entry = iterator.next();
			arr[i++] = new double[]{entry.getKey().getPeriod(), entry.getKey().getIniInventory(), entry.getKey().getIniCash(), entry.getValue()};
		}
		return arr;
	}
	
	
	/**
	 * @param GA
	 * @param GB
	 * @return minus value of resultGB(y*) to resultGA(x)
	 * @date: Jun 2, 2020, 9:52:18 PM 
	 */
	public double[][] getMinusGAGB(double[][] resultGA, double[][] resultGB, double minCash, double fixCost, double variCost){
		int RLength = resultGA.length; 
		int xLength = resultGA[0].length;
		double[][] values = new double[RLength][xLength];
		for (int i = 0; i < RLength; i++) {
			double cash = minCash + i;
			int Qbound = Math.max(0, (int) ((cash - fixCost) / variCost));		
			for (int j = 0; j < xLength; j++) { // this error
				double optimalGB = resultGB[i][j];
				for (int k = j; k < Math.min(j + Qbound + 1, xLength); k++) {
					if (resultGB[i][k] > optimalGB)
						optimalGB = resultGB[i][k];
				}
				values[i][j] = optimalGB - resultGA[i][j] - fixCost;
			}		
		}
		return values;
	}
	
	/**
	 * @param G
	 * @param F
	 * @return minus value of result F to G
	 * @date: Jun 2, 2020, 9:52:18 PM 
	 */
	public double[][] getMinusFG(double[][] resultG, double[][] resultF, double minCash, double fixCost, double variCost){
		int RLength = resultG.length; 
		int xLength = resultG[0].length;
		double[][] values = new double[RLength][xLength];
		for (int i = 0; i < RLength; i++) {
			for (int j = 0; j < xLength; j++) { 
				values[i][j] = resultF[i][j] - resultG[i][j] - variCost * j; // F-G-vx
			}		
		}
		return values;
	}
	
	/**
	 * @param G
	 * @param F
	 * @return values of H
	 * @date: Jun 2, 2020, 9:52:18 PM 
	 */
	public double[][] getH(double[][] resultG, double minCash, double fixCost, double variCost){
		int RLength = resultG.length; 
		int xLength = resultG[0].length;
		double[][] values = new double[RLength][xLength];
		for (int i = 0; i < RLength; i++) {
			double cash = minCash + i;
			for (int j = 0; j < xLength; j++) { 
				int ybound = Math.min(xLength - 1, j + Math.max(0, (int) ((cash - fixCost) / variCost)));	
				double optimalIncre = -fixCost;
				int optyIndex = j;
				for (int k = j + 1; k <= ybound; k++) {
					int cashIndex = (int) (cash - fixCost - variCost * (k - j));
					if (cashIndex < 0)
						continue;
					else {
						if(resultG[cashIndex][k] - resultG[i][j] - fixCost > optimalIncre) { 
							optimalIncre = resultG[cashIndex][k] - resultG[i][j] - fixCost; // be careful
							optyIndex = k;
						}
					}		
				}	
				values[i][j] = optimalIncre; // F-G-vx
			}		
		}
		return values;
	}
	
	
	public boolean checkSingleCrossing(double[][] resultH) {
		int RLength = resultH.length; 
		int xLength = resultH[0].length;
		int lastColumnEnd = 0;
		for (int i = 0; i < RLength; i++) {
			for (int j = xLength - 1; j >= 0; j--) {
				if (resultH[i][j] > -0.1) {
					for (int i1 = i; i1 < RLength; i1++)
						for (int j1 = lastColumnEnd; j1 <= j; j1++) {
							if (resultH[i1][j1] > -0.1)
								continue;
							else 
								return false;
						}
					lastColumnEnd = j + 1;
					break;
				}
			}
		}		
		return true;		
	}
	
	
	/**
	 * @param G
	 * @param F
	 * @return values of H in 3 columns
	 * @date: Jun 2, 2020, 9:52:18 PM 
	 */
	public double[][] getH3Column(double[][] resultG, double minCash, double fixCost, double variCost){
		int RLength = resultG.length; 
		int xLength = resultG[0].length;
		double[][] values = new double[RLength * xLength][3];
		int index = 0;
		for (int i = 0; i < RLength; i++) {
			double cash = minCash + i;
			for (int j = 0; j < xLength; j++) { 
				int ybound = Math.min(xLength - 1, j + Math.max(0, (int) ((cash - fixCost) / variCost)));	
				double optimalIncre = -fixCost;
				int optyIndex = j;
				for (int k = j + 1; k <= ybound; k++) {
					int cashIndex = (int) (cash - fixCost - variCost * (k - j));
					if (cashIndex < 0)
						continue;
					else {
						if(resultG[cashIndex][k] - resultG[i][j] - fixCost > optimalIncre) { 
							optimalIncre = resultG[cashIndex][k] - resultG[i][j] - fixCost; // be careful
							optyIndex = k;
						}
					}		
				}	
				values[index][0] = resultG[i][0];
				values[index][1] = resultG[0][j];
				values[index][2] = optimalIncre; // F-G-vx
				index++;
			}		
		}
		return values;
	}
	
	/**
	 * @param valuse of GB - GA - K
	 * @return check the non-increasing of GB - GA - K for fixed R
	 * @date: Jun 4, 2020, 7:36:55 PM 
	 */
	public boolean checkNonIncreasing(double[][] arr){
		int RLength = arr.length; 
		int xLength = arr[0].length;
		double[][] values = new double[RLength][xLength];
		for (int i = 0; i < RLength; i++) {
			Arrays.fill(values[0], 0);
			for (int j = 0; j < xLength - 1; j++) {
			  if (arr[i][j + 1] - arr[i][j] < 0.1)
				  continue;
			  else {
				  values[i][j] = 1;
				  return false;
			  }
			}
		}		
		return true;
	}
	
	/**
	 * @param valuse of GB - GA - K
	 * @return check the non-decreasing of GB - GA - K for fixed x
	 * @date: Jun 4, 2020, 7:36:55 PM 
	 */
	public boolean checkNonDecreasing(double[][] arr){
		int RLength = arr.length; 
		int xLength = arr[0].length;
		double[][] values = new double[RLength][xLength];
		for (int i = 0; i < RLength; i++) {
			Arrays.fill(values[0], 0);
		}				
		for (int j = 0; j < xLength; j++) {
			for (int i = 0; i < RLength - 1; i++) {
			  if (arr[i + 1][j] - arr[i][j] > -0.1)
				  continue;
			  else {
				  values[i][j] = 1;
				  return false;
			  }
			}
		}		
		return true;
	}
	
	
}
	

