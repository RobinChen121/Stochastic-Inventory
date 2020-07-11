/**
 * @date: Jul 7, 2020
 */
package sdp.cash.multiItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;

import sdp.cash.StateP;
import sdp.cash.StateY;
import sdp.inventory.ImmediateValue.ImmediateValueFunctionV;
import sdp.inventory.StateTransition.StateTransitionFunctionV;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 7, 2020
 * @Desc: Recursion function for V(y1, y2, R) of the two product problem
 *
 */
public class CashRecursionV {
	
	double discountFactor;	
	int TLength;
	GetPmfMulti Pmf;
	double[][][] pmf;
	Map<CashStateMultiYR, double[]> cacheActions = new TreeMap<>();	
	Map<CashStateR, double[]> cacheYStar = new TreeMap<>();	
	Map<CashStateMultiYR, Double> cacheValuesV = new TreeMap<>();
	Map<CashStateMultiYR, Double> cacheValuesPai = new TreeMap<>();
	
	Function<CashStateMultiYR, ArrayList<double[]>> buildActionListV;
	Function<CashStateMultiYR, ArrayList<double[]>> buildActionListPai;
	StateTransitionFunctionV<CashStateMultiYR, double[], CashStateMultiYR> stateTransition;
	ImmediateValueFunctionV<CashStateMultiYR, double[],  Double> immediateValue;
	
	public CashRecursionV(double discountFactor, GetPmfMulti Pmf, Function<CashStateMultiYR, ArrayList<double[]>> buildActionListV,
			Function<CashStateMultiYR, ArrayList<double[]>> buildActionListPai, StateTransitionFunctionV<CashStateMultiYR, double[], CashStateMultiYR> stateTransition, 
			ImmediateValueFunctionV<CashStateMultiYR, double[], Double> immediateValue, int TLength) {
		this.discountFactor = discountFactor;
		this.Pmf = Pmf;
		this.buildActionListV = buildActionListV;
		this.buildActionListPai = buildActionListPai;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		this.TLength = TLength;
		
		// sorted map for recorded actions and values 
		Comparator<CashStateMultiYR> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory1() > o2.getIniInventory1() ? 1 : 
				o1.getIniInventory1() == o2.getIniInventory1() ? o1.getIniInventory2() > o2.getIniInventory2() ? 1 :
					o1.getIniInventory2() == o2.getIniInventory2() ? o1.iniR> o2.iniR  ? 1 :
					o1.iniR == o2.iniR ? 0 : -1: -1 : -1 : -1;
		
		Comparator<CashStateR> keyComparator2 = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
				o1.getPeriod() == o2.getPeriod() ?  o1.iniR> o2.iniR  ? 1 :
							o1.iniR == o2.iniR ? 0 : -1:  -1;
				
		this.cacheActions = new TreeMap<>(keyComparator);
		this.cacheValuesV = new TreeMap<>(keyComparator);		
		this.cacheValuesPai = new TreeMap<>(keyComparator);	
		this.cacheYStar = new TreeMap<>(keyComparator2);
	}
	
	
	
	
	public double getExpectedValuePai(CashStateMultiYR initialState) {
		return this.cacheValuesPai.computeIfAbsent(initialState, s -> {
			int n = s.getPeriod();
			double[][] dAndP = pmf[n - 1]; // demandAndPossibility
			
			double expectValue = 0;
			for (int j = 0; j < dAndP.length; j++) {
				double[] thisDemands = new double[] {dAndP[j][0],  dAndP[j][1]};
				CashStateMultiYR newState = stateTransition.apply(s, thisDemands);
				expectValue += dAndP[j][2] * immediateValue.apply(s, thisDemands);
			}	
			return expectValue;
		});
	}
	
	public double getExpectedValueV(CashStateMultiYR initialState) {
		return this.cacheValuesV.computeIfAbsent(initialState, s -> {
			ArrayList<double[]> actions = buildActionListV.apply(s);
			double val = -Double.MAX_VALUE;
			double[] bestYs = new double[] {0, 0};
			for (int i = 0; i < actions.size(); i++) {
				double[] thisActions = actions.get(i);
				CashStateMultiYR thisState = new CashStateMultiYR(s.getPeriod(), thisActions[0], thisActions[1], s.iniR);
				double thisActionsValue = getExpectedValuePai(thisState);
				
				if (thisActionsValue > val + 0.1) {
					val = thisActionsValue;
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
	* @return optimal y1, y2  for state (x1, x2, R)
	*/
	public double[] getAction(CashStateMultiYR state) {
		return cacheActions.get(state);
	}

	/**
	 * @param period
	 * @param R
	 * @return y1* and y2* for a fixed R in any period
	 * @date: Jul 10, 2020, 12:13:50 PM 
	 */
	public double[] getYStar(CashStateR initialState) {
		return this.cacheYStar.computeIfAbsent(initialState, s -> {
		CashStateMultiYR state = new CashStateMultiYR(s.getPeriod(), 0, 0, s.iniR);
		ArrayList<double[]> actions = buildActionListPai.apply(state);
		double val = -Double.MAX_VALUE;
		double[] bestYs = new double[] {0, 0};
		for (int i = 0; i < actions.size(); i++) {
			double[] thisActions = actions.get(i);
			CashStateMultiYR thisState = new CashStateMultiYR(s.getPeriod(), thisActions[0], thisActions[1], s.iniR);
			double thisActionsValue = getExpectedValuePai(thisState);
			if (thisActionsValue > val + 0.1) {
				val = thisActionsValue;
				bestYs = thisActions;
			}
		}
		return bestYs;			
		});
	}
	
	
	
}
