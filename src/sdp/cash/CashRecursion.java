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
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] QValues = new double[feasibleActions.length];
			double val = optDirection == OptDirection.MIN ? Double.MAX_VALUE
														  : -Double.MAX_VALUE;
		
			double bestOrderQty = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderQty = feasibleActions[i];
				
//				if (s.getPeriod() == 1) { // for debugging
//					orderQty = 110;
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
			arr[i++] =new double[]{entry.getKey().getPeriod(), entry.getKey().getIniInventory(), entry.getKey().getIniCash(), entry.getValue()};
		}
		return arr;
	}
}
