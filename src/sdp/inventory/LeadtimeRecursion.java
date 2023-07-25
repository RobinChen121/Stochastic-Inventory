package sdp.inventory;

import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;

/**
*@author: zhenchen
*@date: Jul 24, 2023, 4:33:44 PM
*@desp: TODO
*
*/

public class LeadtimeRecursion {
	Map<LeadtimeState, Double> cacheActions = new ConcurrentSkipListMap<>();	
	Map<LeadtimeState, Double> cacheValues = new ConcurrentSkipListMap<>();
	
	double[][][] pmf;	
	Function<LeadtimeState, double[]> getFeasibleActions;
	StateTransitionFunction<LeadtimeState, Double, Double, LeadtimeState> stateTransition;
	ImmediateValueFunction<LeadtimeState, Double, Double, Double> immediateValue;
	
	public LeadtimeRecursion(double[][][] pmf, 
	         Function<LeadtimeState, double[]> getFeasibleAction,
	         StateTransitionFunction<LeadtimeState, Double, Double, LeadtimeState> stateTransition,
	         ImmediateValueFunction<LeadtimeState, Double, Double, Double> immediateValue) {
			this.pmf = pmf;
			this.getFeasibleActions = getFeasibleAction;
			this.stateTransition = stateTransition ;
			this.immediateValue = immediateValue;
			
			Comparator<LeadtimeState> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
				o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
					o1.getIniInventory() == o2.getIniInventory() ? o1.getPreQ() > o2.getPreQ() ? 1 :
							o1.getPreQ() == o2.getPreQ() ? 0 : -1 : -1 : -1;
			
			// should have a comparator, or else wrong output
			this.cacheActions = new ConcurrentSkipListMap<>(keyComparator);
			this.cacheValues = new ConcurrentSkipListMap<>(keyComparator);
		}
	
	public double getExpectedValue(LeadtimeState state) {
		return this.cacheValues.computeIfAbsent(state, s -> {						
			double[] feasibleActions = getFeasibleActions.apply(state);
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] QValues = new double[feasibleActions.length];
			double val = Double.MAX_VALUE;

			double bestOrderQty = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderQty = feasibleActions[i]; //
				double thisQValue = 0;								
				for (int j = 0; j < dAndP.length; j++) {
					thisQValue += dAndP[j][1] * immediateValue.apply(s, orderQty, dAndP[j][0]);
					if (s.getPeriod() < pmf.length) {
						LeadtimeState newState = stateTransition.apply(s, orderQty, dAndP[j][0]);
						thisQValue += dAndP[j][1] * getExpectedValue(newState);
					}
				}
				QValues[i] = thisQValue;
				if (QValues[i] < val) {
					val = QValues[i];
					bestOrderQty = orderQty;
				}
			}
			
			this.cacheActions.putIfAbsent(s, bestOrderQty);
			return val;
		});
	}
	
	public double getAction(LeadtimeState LeadtimeState) {
		return cacheActions.get(LeadtimeState);
	}
	
	public Map<LeadtimeState, Double> getCacheActions() {
		return cacheActions;
	}
	
	public Map<LeadtimeState, Double> getCacheValues() {
		return cacheValues;
	}
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(){
		Iterator<Map.Entry<LeadtimeState, Double>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][3];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<LeadtimeState, Double> entry = iterator.next();
			arr[i++] =new double[]{entry.getKey().getPeriod(), entry.getKey().getIniInventory(), entry.getKey().getPreQ(), entry.getValue()};
		}
		return arr;
	}

}


