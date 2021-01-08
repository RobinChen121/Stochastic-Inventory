/**
 * @date: Jul 7, 2020
 */
package sdp.cash.multiItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.function.Function;

import sdp.cash.StateP;
import sdp.cash.StateY;
import sdp.inventory.FinalCash.BoundaryFuncton;
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
	int T;
	GetPmfMulti Pmf;
	double[][][] pmf;
	double[] variCost; 
	Map<CashStateMulti, double[]> cacheActions = new TreeMap<>();	
	Map<CashStateR, double[]> cacheYStar = new TreeMap<>();	
	Map<CashStateR, Double> cacheAlpha = new TreeMap<>();	
	Map<CashStateMulti, Double> cacheValuesV = new TreeMap<>();
	Map<CashStateMultiYR, Double> cacheValuesPai = new TreeMap<>();
	
	Function<CashStateMulti, ArrayList<double[]>> buildActionListV;
	Function<CashStateMultiYR, ArrayList<double[]>> buildActionListPai;
	StateTransitionFunctionV<CashStateMultiYR, double[], CashStateMulti> stateTransition;
	BoundaryFuncton<CashStateMulti, Double> boundFinalCash;
	
	public CashRecursionV(double discountFactor, GetPmfMulti Pmf, Function<CashStateMulti, ArrayList<double[]>> buildActionListV,
			Function<CashStateMultiYR, ArrayList<double[]>> buildActionListPai, StateTransitionFunctionV<CashStateMultiYR, double[], CashStateMulti> stateTransition, 
			BoundaryFuncton<CashStateMulti, Double> boundFinalCash, int T, double[] variCost) {
		this.discountFactor = discountFactor;
		this.Pmf = Pmf;
		this.buildActionListV = buildActionListV;
		this.buildActionListPai = buildActionListPai;
		this.stateTransition = stateTransition;
		this.boundFinalCash = boundFinalCash;
		this.T = T;
		this.variCost = variCost;
		
		// sorted map for recorded yHeads and values 
		Comparator<CashStateMultiYR> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory1() > o2.getIniInventory1() ? 1 : 
				o1.getIniInventory1() == o2.getIniInventory1() ? o1.getIniInventory2() > o2.getIniInventory2() ? 1 :
					o1.getIniInventory2() == o2.getIniInventory2() ? o1.iniR> o2.iniR  ? 1 :
					o1.iniR == o2.iniR ? 0 : -1: -1 : -1 : -1;
					
		Comparator<CashStateMulti> keyComparator1 = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
						o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory1() > o2.getIniInventory1() ? 1 : 
							o1.getIniInventory1() == o2.getIniInventory1() ? o1.getIniInventory2() > o2.getIniInventory2() ? 1 :
								o1.getIniInventory2() == o2.getIniInventory2() ? o1.iniCash > o2.iniCash  ? 1 :
									o1.iniCash == o2.iniCash ? 0 : -1: -1 : -1 : -1;
								
		Comparator<CashStateR> keyComparator2 = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
				o1.getPeriod() == o2.getPeriod() ?  o1.iniR> o2.iniR  ? 1 :
							o1.iniR == o2.iniR ? 0 : -1:  -1;
				
		this.cacheActions = new TreeMap<>(keyComparator1);
		this.cacheValuesV = new TreeMap<>(keyComparator1);		
		this.cacheValuesPai = new TreeMap<>(keyComparator);	
		this.cacheYStar = new TreeMap<>(keyComparator2);
		this.cacheAlpha = new TreeMap<>(keyComparator2);
	}
	
	
	
	
	public double getExpectedValuePai(CashStateMultiYR initialState) {
		return this.cacheValuesPai.computeIfAbsent(initialState, s -> {
			int n = s.getPeriod();
			double[][] dAndP = Pmf.getPmf(n - 1); // demandAndPossibility
			double expectValue = 0;
			for (int j = 0; j < dAndP.length; j++) {
				double[] thisDemands = new double[] {dAndP[j][0],  dAndP[j][1]};
				CashStateMulti newState = stateTransition.apply(s, thisDemands);
				//double thisProfit = immediateValue.apply(newState, thisDemands);
				double thisDemandValue = getExpectedValueV(newState);
				expectValue += dAndP[j][2] * thisDemandValue;
			}	
			return expectValue;
		});
	}
	
	// maybe only stateTransitionFunction needed
	public double getExpectedValueV(CashStateMulti initialState) {
		return this.cacheValuesV.computeIfAbsent(initialState, s -> {
			ArrayList<double[]> yHeads = buildActionListV.apply(s);
			double val = -Double.MAX_VALUE;
			double[] bestYs = new double[] {initialState.getIniInventory1(), initialState.getIniInventory2()};
			if (initialState.getPeriod() <= T) {
				for (int i = 0; i < yHeads.size(); i++) {
					double[] thisActions = yHeads.get(i);
					double iniR = s.iniCash + variCost[0] * s.iniInventory1 + variCost[1] * s.iniInventory2;
					CashStateMultiYR thisState = new CashStateMultiYR(s.getPeriod(), thisActions[0], thisActions[1], iniR);
					double thisActionsValue = getExpectedValuePai(thisState);
//					CashStateR  stateR = new CashStateR(s.getPeriod(), iniR);
//					getYStar(stateR);

					if (thisActionsValue > val + 0.1) {
						val = thisActionsValue;
						bestYs = thisActions;
					}
				}
				double iniR = s.iniCash + variCost[0] * s.iniInventory1 + variCost[1] * s.iniInventory2;
				CashStateR  stateR = new CashStateR(s.getPeriod(), iniR);
				getYStar(stateR);
				this.cacheActions.putIfAbsent(s, bestYs);
			}
			else {
				double finalValue = boundFinalCash.apply(initialState);		
				return finalValue;
			}	
			return val;
		});
	}
	
	/**
	* @Description: return the optimal yHeads for a given state
	* @param @param state
	* @param @return    
	* @return optimal y1, y2  for state (x1, x2, R)
	*/
	public double[] getAction(CashStateMulti state) {
		return cacheActions.get(state);
	}

	/**
	 * @param period
	 * @param R
	 * @return y1* and y2* for a fixed R in any period
	 * @date: Jul 10, 2020, 12:13:50 PM 
	 */
	public double[] getYStar(CashStateR initialState) { // revise
		return this.cacheYStar.computeIfAbsent(initialState, s -> {
		CashStateMultiYR state = new CashStateMultiYR(s.getPeriod(), 0, 0, s.iniR);
		ArrayList<double[]> yHeads = buildActionListPai.apply(state);
		double val = -Double.MAX_VALUE;
		double[] bestYs = new double[] {0, 0};
		for (int i = 0; i < yHeads.size(); i++) {
			double[] thisActions = yHeads.get(i);
			CashStateMultiYR thisState = new CashStateMultiYR(s.getPeriod(), thisActions[0], thisActions[1], s.iniR);
			double thisActionsValue = getExpectedValuePai(thisState);
			
			if (thisActionsValue > val + 0.2) {
				val = thisActionsValue;
				bestYs = thisActions;
			}
//			if ((int)thisActions[0] == 11 && (int)thisActions[1] == 7 && s.period == 1)
//				System.out.println(thisActionsValue);
			
			
		}
//		if (s.period == 1)
//			System.out.println(val);
		if (variCost[0] * bestYs[0] + variCost[1] * bestYs[1] >= s.iniR + 0.1) {
			getAlpha(s); // revise to save computation time
		}
			
		return bestYs;			
		});
	}
	
	// revise
	public double getAlpha(CashStateR initialState) {
		return this.cacheAlpha.computeIfAbsent(initialState, s -> {
			double bestAlpha = 0;
			double bestValue = -Double.MAX_VALUE;
			for (double alpha = 0; alpha <= 1; alpha = alpha + 0.01) {  // stepsize of alpha
				if (initialState.period == 1 && alpha > 0.89)
					System.out.print("");
				double y1 = alpha * s.iniR / variCost[0];
				double y2 = (1 - alpha) * s.iniR / variCost[1];
				CashStateMultiYR state = new CashStateMultiYR(s.period, y1, y2, s.iniR);
				double expectValue = getExpectedValuePai(state);
				if (expectValue > bestValue - 0.1) {
					bestValue = expectValue;
					bestAlpha = alpha;
				}
			}
		return bestAlpha;
		});
	}
	
	/**
	 * 
	 * @return optimal decision table of function Pai
	 */
	public double[][] getOptTable(){
		Iterator<Map.Entry<CashStateR, double[]>> iterator = cacheYStar.entrySet().iterator();
		double[][] arr = new double[cacheYStar.size()][4];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<CashStateR, double[]> entry = iterator.next();
			int period = entry.getKey().getPeriod();
			double R = entry.getKey().getIniR();
			double y1 = entry.getValue()[0];
			double y2 = entry.getValue()[1];
			double c1 = variCost[0]; double c2 = variCost[1];
			double cashConstrained = 0;
			double alpha = 10000;
			if (c1 * y1 + c2 * y2 >= R + 0.1) {
				cashConstrained = 1;
				CashStateR stateR = new CashStateR(period, R);
				alpha = cacheAlpha.get(stateR);
			}
			arr[i++] = new double[]{period, R, variCost[0], variCost[1], y1, y2, cashConstrained, alpha};
		}
		return arr;
	}
	
	
	/**
	 * 
	 * @return detailed optimal decision table including t, x1, x2, w, c1, c2, R, y1, y2, y1*, y2*, cashconstrained, alpha
	 */
	public double[][] getOptTableDetail(){
		Iterator<Map.Entry<CashStateMulti, double[]>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][13]; // revise
		int i = 0;
		Map.Entry<CashStateMulti, double[]> entry;
		while (iterator.hasNext()) {
			entry = iterator.next();
			int period = entry.getKey().period;
			double R = entry.getKey().getIniR(variCost);
			CashStateR stateR = new CashStateR(period, R);
			double[] yStars = cacheYStar.get(stateR);
			double y1 = entry.getValue()[0]; double y2 = entry.getValue()[1];
			
			double cashConstrained = 0;
			double alpha = 10000;
			double x1 = entry.getKey().iniInventory1; double x2 = entry.getKey().iniInventory2;
			double w = entry.getKey().getIniCash(); 
			double[] yHeads = new double[] {0, 0};
			if (x1 < yStars[0]+0.1 && x2 < yStars[1]
					&& variCost[0] * yStars[0] + variCost[1] * yStars[1] < stateR.iniR+0.1) {
				yHeads[0] = yStars[0]; yHeads[1] = yStars[1];
				cashConstrained = 1;
			}
			else if (x1 > yStars[0]-0.1 && x2 > yStars[1]-0.1) {
				yHeads[0] = x1; yHeads[1] = x2;
				cashConstrained = 5;
			}
			else if (x1 > yStars[0]-0.1 && x2 < yStars[1]+0.1) {
				yHeads[0] = x1; yHeads[1] = Math.min(yStars[1], (stateR.iniR - x1*variCost[0]) / variCost[1]);
				cashConstrained = 4;
			}
			else if (x1 < yStars[0]+0.1 && x2 > yStars[1]-0.1){
				yHeads[0] = Math.min(yStars[0], (stateR.iniR -  x2*variCost[1]) / variCost[0]); yHeads[1] = x2;
				cashConstrained = 3;
			}
			else if (x1 < yStars[0]+0.1 && x2 < yStars[1]+0.1
					&& variCost[0] * yStars[0] + variCost[1] * yStars[1] > stateR.iniR-0.1) {
				alpha = cacheAlpha.get(stateR);
				yHeads[0] = alpha * stateR.iniR / variCost[0]; 
				yHeads[1] = (1 - alpha) * stateR.iniR / variCost[1];	
				cashConstrained = 2;
			}				

			arr[i++] = new double[]{period, x1, x2, w, variCost[0], variCost[1], R, yStars[0], yStars[1], cashConstrained, alpha, yHeads[0], yHeads[1]};					
		}
		return arr;
	}
	
	
	/**
	 * 
	 * @return more detailed optimal decision table
	 */
	public double[][] getOptTableDetail2(double[] mean, double[] variance, double[] price, double[] a1, double[] a2){
		Iterator<Map.Entry<CashStateMulti, double[]>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][13]; // revise
		int i = 0;
		Map.Entry<CashStateMulti, double[]> entry;
		while (iterator.hasNext()) {
			entry = iterator.next();
			int period = entry.getKey().period;
			double R = entry.getKey().getIniR(variCost);
			CashStateR stateR = new CashStateR(period, R);
			double[] yStars = cacheYStar.get(stateR);
			double y1 = entry.getValue()[0]; double y2 = entry.getValue()[1];
			
			
			double cashConstrained = 0;
			double alpha = 10000;
			double x1 = entry.getKey().iniInventory1; double x2 = entry.getKey().iniInventory2;
			double w = entry.getKey().getIniCash(); 
			double[] yHeads = new double[] {0, 0};
			if (x1 < yStars[0]+0.1 && x2 < yStars[1]
					&& variCost[0] * yStars[0] + variCost[1] * yStars[1] < stateR.iniR+0.1) {
				yHeads[0] = yStars[0]; yHeads[1] = yStars[1];
				cashConstrained = 1;
			}
			else if (x1 > yStars[0]-0.1 && x2 > yStars[1]-0.1) {
				yHeads[0] = x1; yHeads[1] = x2;
				cashConstrained = 5;
			}
			else if (x1 > yStars[0]-0.1 && x2 < yStars[1]+0.1) {
				yHeads[0] = x1; yHeads[1] = Math.min(yStars[1], (stateR.iniR - x1*variCost[0]) / variCost[1]);
				cashConstrained = 4;
			}
			else if (x1 < yStars[0]+0.1 && x2 > yStars[1]-0.1){
				yHeads[0] = Math.min(yStars[0], (stateR.iniR -  x2*variCost[1]) / variCost[0]); yHeads[1] = x2;
				cashConstrained = 3;
			}
			else if (x1 < yStars[0]+0.1 && x2 < yStars[1]+0.1
					&& variCost[0] * yStars[0] + variCost[1] * yStars[1] > stateR.iniR-0.1) {
				alpha = cacheAlpha.get(stateR);
				yHeads[0] = alpha * stateR.iniR / variCost[0]; 
				yHeads[1] = (1 - alpha) * stateR.iniR / variCost[1];	
				cashConstrained = 2;
			}				

			arr[i++] = new double[]{mean[0], mean[1], variance[0], variance[1], period, x1, x2, w, price[0], price[1], variCost[0], variCost[1], R, yStars[0], yStars[1], cashConstrained, alpha, yHeads[0], yHeads[1], a1[period-1], a2[period-1]};					
		}
		return arr;
	}

		
	
	
}
