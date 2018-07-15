package sdp.inventory;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 12, 2018---10:30:51 AM
*@description:  a recursion class for stochastic dynamic programming
*/

public class Recursion {
	
	Map<State, Double> cacheActions = new TreeMap<>();	
	Map<State, Double> cacheValues = new HashMap<>();
	
	double[][][] pmf;	
	OptDirection optDirection;	
	Function<State, double[]> getFeasibleActions;
	StateTransitionFunction<State, Double, Double, State> stateTransition;
	ImmediateValueFunction<State, Double, Double, Double> immediateValue;
	
	public enum OptDirection{
		MIN,
		MAX
	}
	
	public Recursion(OptDirection optDirection, double[][][] pmf, 
			         Function<State, double[]> getFeasibleAction,
			         StateTransitionFunction<State, Double, Double, State> stateTransition,
			         ImmediateValueFunction<State, Double, Double, Double> immediateValue) {
		this.optDirection = optDirection;
		this.pmf = pmf;
		this.getFeasibleActions = getFeasibleAction;
		this.stateTransition = stateTransition ;
		this.immediateValue = immediateValue;
	}
		
	public StateTransitionFunction<State, Double, Double, State> getStateTransitionFunction(){
		return stateTransition;
	}
	
	public ImmediateValueFunction<State, Double, Double, Double> getImmediateValueFunction(){
		return immediateValue;
	} 
	
	/**
	 * set a tree map for finding s, S, this method can be neglected since default order is same 
	 */
	public void setTreeMapCacheAction() {
		// cacheActions is a sorted map
		Comparator<State> keyComparator = (o1, o2) -> o1.period > o2.period ? 1 : 
			o1.period == o2.period ? o1.initialInventory > o2.initialInventory ? 1 : 
				(o1.initialInventory == o2.initialInventory ? 0 : -1) : -1;
		cacheActions = new TreeMap<>(keyComparator);
	}
	
	
	public double getExpectedValue(State state) {
		return this.cacheValues.computeIfAbsent(state, s -> {			
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
			
			double[] feasibleActions = getFeasibleActions.apply(state);
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] QValues = new double[feasibleActions.length];
			double val = optDirection == OptDirection.MIN ? Double.MAX_VALUE
														  : Double.MIN_VALUE;
			double bestOrderQty = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderQty = feasibleActions[i];
				double thisQValue = 0;								
				for (int j = 0; j < dAndP.length; j++) {
					thisQValue += dAndP[j][1] * immediateValue.apply(s, orderQty, dAndP[j][0]);
					if (s.getPeriod() < pmf.length) {
						State newState = stateTransition.apply(s, orderQty, dAndP[j][0]);
						thisQValue += dAndP[j][1] * getExpectedValue(newState);
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
			this.cacheActions.putIfAbsent(s, bestOrderQty);
			return val;
		});
	}
	
	public double getAction(State state) {
		return cacheActions.get(state);
	}
	
	public Map<State, Double> getCacheActions() {
		return cacheActions;
	}
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(){
		Iterator<Map.Entry<State, Double>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][3];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<State, Double> entry = iterator.next();
			arr[i++] =new double[]{entry.getKey().getPeriod(), entry.getKey().getIniInventory(), entry.getValue()};
		}
		return arr;
	}

}
