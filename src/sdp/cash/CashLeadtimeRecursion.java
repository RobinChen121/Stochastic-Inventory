package sdp.cash;

import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;

/**
*@author: zhenchen
*@date: Mar 2, 2024, 5:28:59 PM
*@desp: TODO
*
*/

public class CashLeadtimeRecursion {
	Map<CashLeadtimeState, Double> cacheActions = new ConcurrentSkipListMap<>();	
	Map<CashLeadtimeState, Double> cacheValues = new ConcurrentSkipListMap<>();
	
	double[][][] pmf;	
	Function<CashLeadtimeState, double[]> getFeasibleActions;
	StateTransitionFunction<CashLeadtimeState, Double, Double, CashLeadtimeState> stateTransition;
	ImmediateValueFunction<CashLeadtimeState, Double, Double, Double> immediateValue;
	
	public CashLeadtimeRecursion(double[][][] pmf, 
	         Function<CashLeadtimeState, double[]> getFeasibleAction,
	         StateTransitionFunction<CashLeadtimeState, Double, Double, CashLeadtimeState> stateTransition,
	         ImmediateValueFunction<CashLeadtimeState, Double, Double, Double> immediateValue) {
			this.pmf = pmf;
			this.getFeasibleActions = getFeasibleAction;
			this.stateTransition = stateTransition;
			this.immediateValue = immediateValue;
			
			Comparator<CashLeadtimeState> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
				o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
					o1.getIniInventory() == o2.getIniInventory() ? o1.getPreQ() > o2.getPreQ() ? 1 :
						o1.getIniCash() == o2.getIniCash() ? o1.getIniCash() > o2.getIniCash() ? 1 :
							o1.getPreQ() == o2.getPreQ() ? 0 : -1 : -1 : -1 : -1;
			
			// should have a comparator, or else wrong output
			this.cacheActions = new ConcurrentSkipListMap<>(keyComparator);
			this.cacheValues = new ConcurrentSkipListMap<>(keyComparator);
		}
	
	public double getExpectedValue(CashLeadtimeState state) {
		return this.cacheValues.computeIfAbsent(state, s -> {						
			double[] feasibleActions = getFeasibleActions.apply(state);
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] QValues = new double[feasibleActions.length];
			double val = -Double.MAX_VALUE;			

			double bestOrderQty = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderQty = feasibleActions[i]; //
//				if (s.getPeriod() == 1) { // only for debugging
//					orderQty = 36;
//				}
				double thisQValue = 0;								
				for (int j = 0; j < dAndP.length; j++) {
					thisQValue += dAndP[j][1] * immediateValue.apply(s, orderQty, dAndP[j][0]);
					if (s.getPeriod() < pmf.length) {
						CashLeadtimeState newState = stateTransition.apply(s, orderQty, dAndP[j][0]);
						thisQValue += dAndP[j][1] * getExpectedValue(newState);
					}
				}
				QValues[i] = thisQValue;
				if (QValues[i] > val) {
					val = QValues[i];
					bestOrderQty = orderQty;
				}
			}
			
			this.cacheActions.putIfAbsent(s, bestOrderQty);
			return val;
		});
	}
	
	public double getAction(CashLeadtimeState LeadtimeState) {
		return cacheActions.get(LeadtimeState);
	}
	
	public Map<CashLeadtimeState, Double> getCacheActions() {
		return cacheActions;
	}
	
	public Map<CashLeadtimeState, Double> getCacheValues() {
		return cacheValues;
	}
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(){
		Iterator<Map.Entry<CashLeadtimeState, Double>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][3];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<CashLeadtimeState, Double> entry = iterator.next();
			arr[i++] =new double[]{entry.getKey().getPeriod(), entry.getKey().getIniInventory(), entry.getKey().getIniCash(), entry.getKey().getPreQ(), entry.getValue()};
		}
		return arr;
	}
	
	
	
	
	
	
	

}


