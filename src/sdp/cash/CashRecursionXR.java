package sdp.cash;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Feb 22, 2020---14:40:05 PM
*@description:  recursion class for cash flow problem with state x and R
*/

public class CashRecursionXR {
	Map<CashStateXR, Double> cacheActions = new  ConcurrentSkipListMap<>();	
	Map<CashStateXR, Double> cacheValues = new  ConcurrentSkipListMap<>();
	
	double[][][] pmf;	
	OptDirection optDirection;	
	Function<CashStateXR, double[]> getFeasibleActions;
	StateTransitionFunction<CashStateXR, Double, Double, CashStateXR> stateTransition;
	ImmediateValueFunction<CashStateXR, Double, Double, Double> immediateValue;
	double discountFactor;
	
	public enum OptDirection{
		MIN,
		MAX
	}
	
	public CashRecursionXR(OptDirection optDirection, double[][][] pmf, 
			         Function<CashStateXR, double[]> getFeasibleAction,
			         StateTransitionFunction<CashStateXR, Double, Double, CashStateXR> stateTransition,
			         ImmediateValueFunction<CashStateXR, Double, Double, Double> immediateValue, 
			         double discountFactor) {
		this.optDirection = optDirection;
		this.pmf = pmf;
		this.getFeasibleActions = getFeasibleAction;
		this.stateTransition = stateTransition ;
		this.immediateValue = immediateValue;
		Comparator<CashStateXR> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? o1.iniR> o2.iniR  ? 1 :
					o1.iniR == o2.iniR ? 0 : -1 : -1 : -1;
		this.cacheActions = new ConcurrentSkipListMap<>(keyComparator);
		this.cacheValues = new ConcurrentSkipListMap<>(keyComparator);
		this.discountFactor = discountFactor;
	}
		
	public StateTransitionFunction<CashStateXR, Double, Double, CashStateXR> getStateTransitionFunction(){
		return stateTransition;
	}
	
	public ImmediateValueFunction<CashStateXR, Double, Double, Double> getImmediateValueFunction(){
		return immediateValue;
	} 
	
	/**
	 * set a tree map 
	 */
	public void setTreeMapCacheAction() {
		// cacheActions is a sorted map
		Comparator<CashStateXR> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? o1.iniR > o2.iniR  ? 1 :
					o1.iniR == o2.iniR ? 0 : -1 : -1 : -1;
		cacheActions = new TreeMap<>(keyComparator);
	}
	
	
	public double getExpectedValue(CashStateXR initialState) {
		return this.cacheValues.computeIfAbsent(initialState, s -> {			
					
			double[] feasibleActions = getFeasibleActions.apply(initialState);
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] YValues = new double[feasibleActions.length];
			double val = optDirection == OptDirection.MIN ? Double.MAX_VALUE
														  : -Double.MAX_VALUE;
		
			double bestY = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderY = feasibleActions[i];				
				double thisYValue = 0;								
				for (int j = 0; j < dAndP.length; j++) {
					//System.out.println(dAndP[j][0]);
					double thisValue = immediateValue.apply(s, orderY, dAndP[j][0]);
					thisYValue += dAndP[j][1] * thisValue;
					if (s.getPeriod() < pmf.length) {
						double a = dAndP[j][0];
						CashStateXR newState = stateTransition.apply(s, orderY, dAndP[j][0]);
						if (dAndP[j][0] < 0)
							System.out.println(a);
						thisYValue += dAndP[j][1] * discountFactor * getExpectedValue(newState);
					}
				}
				YValues[i] = thisYValue;
				if (optDirection == OptDirection.MIN) {
					if (YValues[i] < val) {
						val = YValues[i];
						bestY = orderY;
					}
				}
				else {
					if (YValues[i] > val) {
						val = YValues[i];
						bestY = orderY;
					}
				}
			}
			try {
			this.cacheActions.putIfAbsent(s, bestY);
			}
			catch (Exception e) {
				System.out.println("error");
			}
			return val;
		});
	}
	
	public double getAction(CashStateXR state) {
		return cacheActions.get(state);
	}
	
	public Map<CashStateXR, Double> getCacheActions() {
		return cacheActions;
	}
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(){
		Iterator<Map.Entry<CashStateXR, Double>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][3];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<CashStateXR, Double> entry = iterator.next();
			double a = entry.getKey().getIniInventory();
			double b = entry.getKey().getIniR();
			double S = b - entry.getKey().unitVariCost * a;
			arr[i++] =new double[]{entry.getKey().getPeriod(), entry.getKey().getIniInventory(), S, entry.getKey().getIniR(), entry.getValue()};
		}
		return arr;
	}
}
