/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 24, 2019, 5:59:55 PM
 * @Desc: recursion function for multi item for states: x1, x2, R
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
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;





public class CashRecursionMultiXR {	
	double discountFactor;	
	int TLength;
	GetPmfMulti Pmf;
	double[][][] pmf;
	Map<CashStateMultiXR, double[]> cacheActions = new ConcurrentSkipListMap<>();	
	Map<CashStateMultiXR, Double> cacheValues = new ConcurrentSkipListMap<>();
	
	Function<CashStateMultiXR, ArrayList<double[]>> buildActionList;
	StateTransitionFunction<CashStateMultiXR, double[], double[], CashStateMultiXR> stateTransition;
	ImmediateValueFunction<CashStateMultiXR, double[], double[], Double> immediateValue;
	
	public CashRecursionMultiXR(double discountFactor, GetPmfMulti Pmf, Function<CashStateMultiXR, ArrayList<double[]>> buildActionList,
			StateTransitionFunction<CashStateMultiXR, double[], double[], CashStateMultiXR> stateTransition, 
			ImmediateValueFunction<CashStateMultiXR, double[], double[], Double> immediateValue, int TLength) {
		this.discountFactor = discountFactor;
		this.Pmf = Pmf;
		this.buildActionList = buildActionList;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		this.TLength = TLength;
		
		// sorted map for recorded actions and values 
		Comparator<CashStateMultiXR> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory1() > o2.getIniInventory1() ? 1 : 
				o1.getIniInventory1() == o2.getIniInventory1() ? o1.getIniInventory2() > o2.getIniInventory2() ? 1 :
					o1.getIniInventory2() == o2.getIniInventory2() ? o1.iniR> o2.iniR  ? 1 :
					o1.iniR == o2.iniR ? 0 : -1: -1 : -1 : -1;
		this.cacheActions = new ConcurrentSkipListMap<>(keyComparator);
		this.cacheValues = new ConcurrentSkipListMap<>(keyComparator);		
	}
	
	

	public double getExpectedValue(CashStateMultiXR initialState) {
		return this.cacheValues.computeIfAbsent(initialState, s -> {
			ArrayList<double[]> actions = buildActionList.apply(s);
			double[][] dAndP = Pmf.getPmf(s.getPeriod() - 1);
			double val = -Double.MAX_VALUE;
			
//			double pSum = 0;
//			for (int i = 0; i < dAndP.length; i++) {
//				pSum += dAndP[i][2];
//			}
//			System.out.println(pSum);
			
			double[] actionValues = new double[actions.size()];
			double[] bestYs = new double[] {0, 0};
			for (int i = 0; i < actions.size(); i++) {
				double[] thisActions = actions.get(i);
				double thisActionsValue = 0;
				for (int j = 0; j < dAndP.length; j++) {
					double[] thisDemands = new double[] {dAndP[j][0],  dAndP[j][1]};
					thisActionsValue += dAndP[j][2] * immediateValue.apply(s, thisActions, thisDemands);
					if (s.getPeriod()  < TLength) {
						CashStateMultiXR newState = stateTransition.apply(s, thisActions, thisDemands);
						thisActionsValue += dAndP[j][2] * discountFactor * getExpectedValue(newState);
					}
				}
				actionValues[i] = thisActionsValue;
				if (actionValues[i] > val + 0.1) {
					val = actionValues[i];
					bestYs = thisActions;
				}
			}
			this.cacheActions.putIfAbsent(s, bestYs);
			return val;
		});
	}
	
	
	
	/**
	* @Description: return the optimal actions for a given state
	* @param @param state
	* @param @return    
	* @return y1, y2  
	*/
	public double[] getAction(CashStateMultiXR state) {
		return cacheActions.get(state);
	}
	
	
	/**
	* @Description: return total optimal actions for all states
	* @param @return    
	* @return Map<CashStateMulti,Actions>   
	*/
	public Map<CashStateMultiXR, double[]> getCacheActions() {
		return cacheActions;
	}
	
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(double[] variCost){
		Iterator<Map.Entry<CashStateMultiXR, double[]>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][10];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<CashStateMultiXR, double[]> entry = iterator.next();
			double period = entry.getKey().getPeriod();
			double x1 = entry.getKey().getIniInventory1();
			double x2 = entry.getKey().getIniInventory2();
			double R = entry.getKey().getIniR();
			double y1 = entry.getValue()[0];
			double y2 = entry.getValue()[1];
			double Q1 = y1 - x1;
			double Q2 = y2 - x2;
			double alpha = 10000;
			double boolAlpha = 0;
			double w = R - x1 * variCost[0] - x2 * variCost[1];
			if (R <= variCost[0] * y1 + variCost[1] * y2 + 0.1 && Q1 > 0.1 && Q2 > 0.1) {
				boolAlpha = 1;
				alpha = variCost[0] * y1 / R;
			}			
			arr[i++] = new double[]{period, x1, x2, w, R, boolAlpha, alpha, y1, y2, variCost[0], variCost[1]};
		}
		return arr;
	}

}
