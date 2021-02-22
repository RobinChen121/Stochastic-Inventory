/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jan 20, 2020, 5:19:20 PM
 * @Desc:  Recursion function to compute G
 *
 *
 * 
 */
package sdp.cash.multiItem;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;

import sdp.inventory.State;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;


public class RecursionG {
	Map<State, Double> cacheActions = new TreeMap<>();	
	Map<State, Double> cacheValues = new TreeMap<>();
	Map<State, double[]> cacheActionValues = new TreeMap<>();
	
	double[][][] pmf;		
	Function<State, double[]> getFeasibleActions;
	StateTransitionFunction<State, Double, Double, State> stateTransition;
	ImmediateValueFunction<State, Double, Double, Double> immediateValue;
	
	double tOptY[]; // optimal Y in each period
	
	public RecursionG( double[][][] pmf, Function<State, double[]> getFeasibleAction,
			StateTransitionFunction<State, Double, Double, State> stateTransition,
			ImmediateValueFunction<State, Double, Double, Double> immediateValue) {
		this.pmf = pmf;
		this.getFeasibleActions = getFeasibleAction;
		this.stateTransition = stateTransition ;
		this.immediateValue = immediateValue;
		Comparator<State> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? 0 : -1 : -1;
		this.cacheActions = new TreeMap<>(keyComparator);
		this.cacheValues = new TreeMap<>(keyComparator);
		this.cacheActionValues = new TreeMap<>(keyComparator);
		this.tOptY = new double[pmf.length + 1];
		Arrays.fill(tOptY, -100); // initialize tOptY
	}
	
	public StateTransitionFunction<State, Double, Double, State> getStateTransitionFunction(){
		return stateTransition;
	}
	
	public ImmediateValueFunction<State, Double, Double, Double> getImmediateValueFunction(){
		return immediateValue;
	}
	
	public double getExpectedValue(State state) {
		return this.cacheValues.computeIfAbsent(state, s -> {	
			double[] feasibleActions = getFeasibleActions.apply(state);
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] QValues = new double[feasibleActions.length];
			double val = -Double.MAX_VALUE;
			double bestOrderQty = 0;
			double lastBestY = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderQty = feasibleActions[i];
				double thisQValue = 0;								
				for (int j = 0; j < dAndP.length; j++) {
					double thisValue = immediateValue.apply(s, orderQty, dAndP[j][0]);
					thisQValue += dAndP[j][1] * immediateValue.apply(s, orderQty, dAndP[j][0]);
					if (s.getPeriod() < pmf.length) {
						lastBestY = getLastOptY(state);
						double maxY = Math.max(s.getIniInventory() + orderQty - dAndP[j][0], lastBestY);
						State newState = new State(s.getPeriod() + 1, maxY);						
						thisQValue += dAndP[j][1] * getExpectedValue(newState);
					}
					else {
						tOptY[pmf.length] = 0;
					}
				}
				QValues[i] = thisQValue;
				if (QValues[i] > val) {
					val = QValues[i];
					bestOrderQty = orderQty;
				}
			}
			
			if (s.getPeriod() == 3)
				System.out.println();
			this.cacheActions.putIfAbsent(s, bestOrderQty);
			this.cacheActionValues.putIfAbsent(s, new double[] {bestOrderQty, val});
			return val;
		});
	}
	
	public double getLastOptY(State state) {
		double lastOptY = 0;
		if (state.getPeriod() == pmf.length)
			lastOptY = 0;
		else {
			if (state.getPeriod() == 3)
				System.out.println();
			if (tOptY[state.getPeriod()] < 0) {
				double optValue = -Double.MAX_VALUE;
				Iterator<Map.Entry<State, double[]>> iterator = this.cacheActionValues.entrySet().iterator();
				while (iterator.hasNext()) {
					Map.Entry<State, double[]> entry = iterator.next();
					if (entry.getKey().getPeriod() == state.getPeriod() + 1) {
						if (entry.getValue()[1] > optValue) {
							optValue = entry.getValue()[1];
							lastOptY = entry.getValue()[0] + entry.getKey().getIniInventory();
							tOptY[state.getPeriod()] = lastOptY;
						}
					}
					if (entry.getKey().getPeriod() == state.getPeriod() + 2)
						break;
				}
			}
		}
		
		return lastOptY;
	}
	
	
	public double[] getOptY() {
		return tOptY;
	}
	
	public double getAction(State state) {
		return cacheActions.get(state);
	}
}
