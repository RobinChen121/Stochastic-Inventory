/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 24, 2019, 5:59:55 PM
 * @Desc: 
 *
 *
 * 
 */
package sdp.cash.multiItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;

import cash.multiItem.ImmediateValueFunction;
import cash.multiItem.StateTransitionFunction;



public class CashRecursionMulti {	
	double discountFactor;	
	int TLength;
	GetPmfMulti Pmf;
	Map<CashStateMulti, Actions> cacheActions = new TreeMap<>();	
	Map<CashStateMulti, Double> cacheValues = new TreeMap<>();
	
	Function<CashStateMulti, ArrayList<Actions>> buildActionList;
	StateTransitionFunction<CashStateMulti, Actions, Demands, CashStateMulti> stateTransition;
	ImmediateValueFunction<CashStateMulti, Actions, Demands, Double> immediateValue;
	
	public CashRecursionMulti(double discountFactor, GetPmfMulti Pmf, Function<CashStateMulti, ArrayList<Actions>> buildActionList,
			StateTransitionFunction<CashStateMulti, Actions, Demands, CashStateMulti> stateTransition, 
			ImmediateValueFunction<CashStateMulti, Actions, Demands, Double> immediateValue, int TLength) {
		this.discountFactor = discountFactor;
		this.Pmf = Pmf;
		this.buildActionList = buildActionList;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		this.TLength = TLength;
		
		// sorted map for recorded actions and values 
		Comparator<CashStateMulti> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory1() > o2.getIniInventory1() ? 1 : 
				o1.getIniInventory1() == o2.getIniInventory1() ? o1.getIniInventory2() > o2.getIniInventory2() ? 1 :
					o1.getIniInventory2() == o2.getIniInventory2() ? o1.iniCash > o2.iniCash  ? 1 :
					o1.iniCash == o2.iniCash ? 0 : -1: -1 : -1 : -1;
		this.cacheActions = new TreeMap<>(keyComparator);
		this.cacheValues = new TreeMap<>(keyComparator);		
	}
	
	
	public double getExpectedValue(CashStateMulti initialState) {
		return this.cacheValues.computeIfAbsent(initialState, s -> {
			ArrayList<Actions> actions = buildActionList.apply(s);
			double[][] dAndP = Pmf.getPmf(s.getPeriod() - 1);
			double val = -Double.MAX_VALUE;
			
//			double pSum = 0;
//			for (int i = 0; i < dAndP.length; i++) {
//				pSum += dAndP[i][2];
//			}
//			System.out.println(pSum);
			
			double[] actionValues = new double[actions.size()];
			Actions bestActions = new Actions(0, 0);
			for (int i = 0; i < actions.size(); i++) {
				Actions thisActions = actions.get(i);
				double thisActionsValue = 0;
				for (int j = 0; j < dAndP.length; j++) {
					Demands thisDemands = new Demands((int) dAndP[j][0], (int) dAndP[j][1]);
					thisActionsValue += dAndP[j][2] * immediateValue.apply(s, thisActions, thisDemands);
					if (s.getPeriod()  < TLength) {
						CashStateMulti newState = stateTransition.apply(s, thisActions, thisDemands);
						thisActionsValue += dAndP[j][2] * discountFactor * getExpectedValue(newState);
					}
				}
				actionValues[i] = thisActionsValue;
				if (actionValues[i] > val + 0.1) {
					val = actionValues[i];
					bestActions = new Actions(thisActions.getFirstAction(), thisActions.getSecondAction());
				}
			}
			this.cacheActions.putIfAbsent(s, bestActions);
			return val;
		});
	}
	
	
	/**
	* @Description: return the optimal actions for a given state
	* @param @param state
	* @param @return    
	* @return Actions   
	*/
	public Actions getAction(CashStateMulti state) {
		return cacheActions.get(state);
	}
	
	
	/**
	* @Description: return total optimal actions for all states
	* @param @return    
	* @return Map<CashStateMulti,Actions>   
	*/
	public Map<CashStateMulti, Actions> getCacheActions() {
		return cacheActions;
	}
	
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(double[] variCost){
		Iterator<Map.Entry<CashStateMulti, Actions>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][10];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<CashStateMulti, Actions> entry = iterator.next();
			double period = entry.getKey().getPeriod();
			double x1 = entry.getKey().getIniInventory1();
			double x2 = entry.getKey().getIniInventory2();
			double w = entry.getKey().getIniCash();
			double Q1 = entry.getValue().getFirstAction();
			double Q2 = entry.getValue().getSecondAction();
			double alpha = 10000;
			double boolAlpha = 0;
			double R = w + x1 * variCost[0] + x2 * variCost[1];
			if (w <= variCost[0] * Q1 + variCost[1] * Q2 && Q1 > 0 && Q2 > 0) {
				boolAlpha = 1;
				alpha = variCost[0] * Q1 / R;
			}			
			arr[i++] = new double[]{period, x1, x2, w, R, boolAlpha, alpha, Q1, Q2, variCost[0], variCost[1]};
		}
		return arr;
	}

}
