package sdp.cash;

import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;

import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: June 11, 2023---11:17:05 AM
*@description:  recursion class for survival maximization problem
*/

public class RiskRecursion {
	Map<RiskState, Double> cacheActions = new  ConcurrentSkipListMap<>();	
	Map<RiskState, Double> cacheValues = new  ConcurrentSkipListMap<>();
	
	double[][][] pmf;	
	Function<RiskState, double[]> getFeasibleActions;
	StateTransitionFunction<RiskState, Double, Double, RiskState> stateTransition;
	ImmediateValueFunction<RiskState, Double, Double, Double> immediateValue;
	
	
	public RiskRecursion(double[][][] pmf, 
			         Function<RiskState, double[]> getFeasibleAction,
			         StateTransitionFunction<RiskState, Double, Double, RiskState> stateTransition,
			         ImmediateValueFunction<RiskState, Double, Double, Double> immediateValue) {
		this.pmf = pmf;
		this.getFeasibleActions = getFeasibleAction;
		this.stateTransition = stateTransition ;
		this.immediateValue = immediateValue;
		Comparator<RiskState> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? o1.iniCash > o2.iniCash  ? 1 :
					o1.iniCash == o2.iniCash ? 0 : -1 : -1 : -1;
		this.cacheActions = new  ConcurrentSkipListMap<>(keyComparator);
		this.cacheValues = new  ConcurrentSkipListMap<>(keyComparator);
	}
		
	public StateTransitionFunction<RiskState, Double, Double, RiskState> getStateTransitionFunction(){
		return stateTransition;
	}
	
	public ImmediateValueFunction<RiskState, Double, Double, Double> getImmediateValueFunction(){
		return immediateValue;
	}
	
	public void setTreeMapCacheAction() {
		// cacheActions is a sorted map
		Comparator<RiskState> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? o1.iniCash > o2.iniCash  ? 1 :
					o1.iniCash == o2.iniCash ? 0 : -1 : -1 : -1;
		cacheActions = new TreeMap<>(keyComparator);
	}
	
	public double getSurvProb(RiskState initialState) {
		return this.cacheValues.computeIfAbsent(initialState, s -> {			
			double[] feasibleActions = getFeasibleActions.apply(initialState);
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] QValues = new double[feasibleActions.length];
			double val = -Double.MAX_VALUE;
		
			double bestOrderQty = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderQty = feasibleActions[i];
//				if (s.getPeriod() == 1)
//					orderQty = 8;
				
				double thisQProb = 0;								
				for (int j = 0; j < dAndP.length; j++) {
					double randomDemand = dAndP[j][0];
					double dProb = dAndP[j][1];
					if (s.getPeriod() == pmf.length) {
						double thisDFinalCash = s.iniCash + immediateValue.apply(s, orderQty, randomDemand);
						double thisDProb = thisDFinalCash >= 0 ? 1 : 0;	
						thisQProb += dProb * thisDProb;
					}
					if (s.getPeriod() < pmf.length) { 
						RiskState newState = stateTransition.apply(s, orderQty, dAndP[j][0]);
						double thisDProb = 0;
						if (newState.iniCash < 0) {
							thisDProb = 0;
						}
						else {
							thisDProb = getSurvProb(newState);
						}
						thisQProb += dAndP[j][1] * thisDProb;
					}
				}
				QValues[i] = thisQProb;

				if (QValues[i] > val) {
					val = QValues[i];
					bestOrderQty = orderQty;
				}
			}
			this.cacheActions.putIfAbsent(s, bestOrderQty);
			return val;
		});
	}
	
	
	public double getAction(RiskState state) {
		return cacheActions.get(state);
	}
	
	public Map<RiskState, Double> getCacheActions() {
		return cacheActions;
	}
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(){
		Iterator<Map.Entry<RiskState, Double>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][4];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<RiskState, Double> entry = iterator.next();
			double temp = entry.getKey().getBankruptBefore() ? 1 : 0;
			arr[i++] = new double[]{entry.getKey().getPeriod(), entry.getKey().getIniInventory(), entry.getKey().getIniCash(), temp, entry.getValue()};
		}
		return arr;
	}	
	
}
