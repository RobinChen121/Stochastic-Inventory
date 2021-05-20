package sdp.cash.multiItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;


import sdp.inventory.FinalCash.BoundaryFuncton;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.inventory.StateTransition.StateTransitionFunctionV;

/**
 * @author chen
 * @email: 15011074486@163.com
 * @Date: 2021 Feb 26, 11:27:35
 * @Description: cash recursion for function V(x1, x2, w) for the two product problem
 * 
 */
public class CashRecursionV2 {
	double discountFactor;	
	int T;
	GetPmfMulti Pmf;
	double[][][] pmf;
	double[] variCost; 
	Map<CashStateMulti, double[]> cacheActions = new ConcurrentSkipListMap<>();	
	Map<CashStateMulti, double[]> cacheYStar = new ConcurrentSkipListMap<>();	
	Map<CashStateMulti, Double> cacheAlpha = new ConcurrentSkipListMap<>();	
	Map<CashStateMulti, Double> cacheValuesV = new ConcurrentSkipListMap<>();
	Map<CashStateMultiXYW, Double> cacheValuesPai = new ConcurrentSkipListMap<>();
	
	Function<CashStateMulti, ArrayList<double[]>> buildActionListV;
	Function<CashStateMulti, ArrayList<double[]>> buildActionListPai;
	StateTransitionFunction<CashStateMulti, double[], double[], CashStateMulti> stateTransition;
	BoundaryFuncton<CashStateMulti, Double> boundFinalCash;
	
	public CashRecursionV2(double discountFactor, GetPmfMulti Pmf, Function<CashStateMulti, ArrayList<double[]>> buildActionListV,
			Function<CashStateMulti, ArrayList<double[]>> buildActionListPai, StateTransitionFunction<CashStateMulti, double[], double[], CashStateMulti> stateTransition, 
			BoundaryFuncton<CashStateMulti, Double> boundFinalCash, int T, double[] variCost) {
		this.discountFactor = discountFactor;
		this.Pmf = Pmf;
		this.buildActionListV = buildActionListV;
		this.buildActionListPai = buildActionListPai;
		this.stateTransition = stateTransition;
		this.boundFinalCash = boundFinalCash;
		this.T = T;
		this.variCost = variCost;
		
		Comparator<CashStateMulti> keyComparator1 = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory1() > o2.getIniInventory1() ? 1 : 
				o1.getIniInventory1() == o2.getIniInventory1() ? o1.getIniInventory2() > o2.getIniInventory2() ? 1 :
					o1.getIniInventory2() == o2.getIniInventory2() ? o1.iniCash > o2.iniCash  ? 1 :
						o1.iniCash == o2.iniCash ? 0 : -1: -1 : -1 : -1;
					
		Comparator<CashStateMultiXYW> keyComparator2 = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
				o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory1() > o2.getIniInventory1() ? 1 : 
					o1.getIniInventory1() == o2.getIniInventory1() ? o1.getIniInventory2() > o2.getIniInventory2() ? 1 :
						o1.getIniInventory2() == o2.getIniInventory2() ? o1.getY1() > o2.getY1() ? 1 :
								o1.getY1() == o2.getY1() ? o1.getY2() > o2.getY2() ? 1 :
								o1.getY2() == o2.getY2() ? o1.iniW > o2.iniW ? 1 :
									o1.iniW == o2.iniW ? 0 : -1: -1 : -1 : -1 : -1 : -1;
					
		this.cacheActions = new ConcurrentSkipListMap<>(keyComparator1);
		this.cacheValuesV = new ConcurrentSkipListMap<>(keyComparator1);		
		this.cacheValuesPai = new ConcurrentSkipListMap<>(keyComparator2);	
		this.cacheYStar = new ConcurrentSkipListMap<>(keyComparator1);
		this.cacheAlpha = new ConcurrentSkipListMap<>(keyComparator1);		
	}
	
	
	public double getExpectedValuePai(CashStateMulti initialState, double[] actions) {
		CashStateMultiXYW state = new CashStateMultiXYW(initialState.period, initialState.iniInventory1, initialState.iniInventory2, actions[0], actions[1], initialState.iniCash);
		return this.cacheValuesPai.computeIfAbsent(state, s -> {
			int n = s.getPeriod();
			double[][] dAndP = Pmf.getPmf(n - 1); // demandAndPossibility
			double expectValue = 0;
			for (int j = 0; j < dAndP.length; j++) {
				double[] thisDemands = new double[] {dAndP[j][0],  dAndP[j][1]};
				//double[] actions = new double[] {s.y1, s.y2};
				//CashStateMulti state = new CashStateMulti(s.period, s.iniInventory1, s.iniInventory2, s.iniW);
				CashStateMulti newState = stateTransition.apply(initialState, actions, thisDemands);
				double thisDemandValue = getExpectedValueV(newState);
				expectValue += dAndP[j][2] * thisDemandValue;
			}	
			return expectValue;
		});
	}
	
	public double getExpectedValueV(CashStateMulti initialState) {
		return this.cacheValuesV.computeIfAbsent(initialState, s -> {
			double val = -Double.MAX_VALUE;
			double valYstar = -Double.MAX_VALUE;
			ArrayList<double[]> yHeads = buildActionListV.apply(s);
			ArrayList<double[]> ystars = buildActionListPai.apply(s);
			double[][] dAndP = Pmf.getPmf(s.getPeriod() - 1);
			double[] bestYheads = new double[] { s.iniInventory1, s.iniInventory2 };			
			double[] bestYStars = new double[] { s.iniInventory1, s.iniInventory2 };		

			for (int i = 0; i < ystars.size(); i++) {
				double[] thisActions = ystars.get(i);
//				if (s.getPeriod() == 2 && s.iniCash > 51)
//					thisActions = new double[] {11, 11};
				double thisActionsValue = 0;
				for (int j = 0; j < dAndP.length; j++) {
					double[] thisDemands = new double[] { dAndP[j][0], dAndP[j][1] };
					CashStateMulti newState = stateTransition.apply(s, thisActions, thisDemands);
					if (s.getPeriod() < T)
						thisActionsValue += dAndP[j][2] * discountFactor * getExpectedValueV(newState);
					else
						thisActionsValue += dAndP[j][2] * discountFactor * boundFinalCash.apply(newState);
				}
				if (variCost[0] * thisActions[0] + variCost[1] * thisActions[1] < s.iniCash + 0.1) { // for computing y heads
					if (thisActionsValue > val ) {
						val = thisActionsValue;
						bestYheads = thisActions;
						valYstar = thisActionsValue;
						bestYStars = thisActions;
					}		
				}
				if (variCost[0] * thisActions[0] + variCost[1] * thisActions[1] >= s.iniCash + 0.1) { // for computing y stars
					if (thisActionsValue > valYstar ) {
						valYstar = thisActionsValue;
						bestYStars = thisActions;
					}		
				}
//				if (s.getPeriod() == 2 && s.iniCash > 51)
//					System.out.println(thisActionsValue);
			}		
			
			if (variCost[0] * (bestYStars[0] - s.iniInventory1) + variCost[1] * (bestYStars[1] - s.iniInventory2) >= s.iniCash + 0.1) {
				double bestAlpha = 0;
				double bestValue = -Double.MAX_VALUE;
				for (double alpha = 0; alpha <= 1; alpha = alpha + 0.1) {  // stepsize of alpha
					double y1 = alpha * s.iniCash / variCost[0] + s.iniInventory1;
					double y2 = (1 - alpha) * s.iniCash / variCost[1] + s.iniInventory2;
					double[] thisActions =  new double[] {y1, y2};
					double thisActionsValue = 0;
					for (int j = 0; j < dAndP.length; j++) {
						double[] thisDemands = new double[] { dAndP[j][0], dAndP[j][1] };
						CashStateMulti newState = stateTransition.apply(s, thisActions, thisDemands);
						if (s.getPeriod() < T)
							thisActionsValue += dAndP[j][2] * discountFactor * getExpectedValueV(newState);
						else
							thisActionsValue += dAndP[j][2] * discountFactor * boundFinalCash.apply(newState);
					}
					if (thisActionsValue > bestValue - 0.1) {
						bestValue = thisActionsValue;
						bestAlpha = alpha;
					}
				}
				cacheAlpha.putIfAbsent(initialState, bestAlpha);
			}
			
			this.cacheYStar.putIfAbsent(s, bestYStars);
			this.cacheActions.putIfAbsent(s, bestYheads);	
			return val;
		});
	}

	
	/**
	* @Description: return the optimal yHeads for a given state
	* @param @param state
	* @param @return    
	* @return optimal y1, y2  for state (x1, x2, w)
	*/
	public double[] getAction(CashStateMulti state) {
		return cacheActions.get(state);
	}
	
	public double getAlpha(CashStateMulti state) {
		return cacheAlpha.get(state);
	}
	
	/**
	 * @param period
	 * @param state(x1, x2, w)
	 * @return y1* and y2* for a fixed (x1, x2, w) in any period
	 * @date: Feb 26, 2021, 14:13:50 PM 
	 */
	public double[] getYStar(CashStateMulti initialState) { 
		return this.cacheYStar.computeIfAbsent(initialState, s -> {
			ArrayList<double[]> actions = buildActionListPai.apply(initialState);
			double val = -Double.MAX_VALUE;
			double[][] dAndP = Pmf.getPmf(s.getPeriod() - 1);
			double[] bestYStar = new double[] { s.iniInventory1, s.iniInventory2 };
			for (int i = 0; i < actions.size(); i++) {
				double[] thisActions = actions.get(i);
				
				double thisActionsValue = 0;
				for (int j = 0; j < dAndP.length; j++) {
					double[] thisDemands = new double[] { dAndP[j][0], dAndP[j][1] };
					CashStateMulti newState = stateTransition.apply(s, thisActions, thisDemands);
					if (s.getPeriod() < T)
						thisActionsValue += dAndP[j][2] * discountFactor * getExpectedValueV(newState);
					else
						thisActionsValue += dAndP[j][2] * discountFactor * boundFinalCash.apply(newState);
				}
				if (thisActionsValue > val + 0.1) {
					val = thisActionsValue;
					bestYStar = thisActions;
				}

			}

			if (variCost[0] * (bestYStar[0] - s.iniInventory1) + variCost[1] * (bestYStar[1] - s.iniInventory2) >= s.iniCash + 0.1) {
				double bestAlpha = 0;
				double bestValue = -Double.MAX_VALUE;
				for (double alpha = 0; alpha <= 1; alpha = alpha + 0.01) {  // stepsize of alpha
					double y1 = alpha * s.iniCash / variCost[0] + s.iniInventory1;
					double y2 = (1 - alpha) * s.iniCash / variCost[1] + s.iniInventory2;
					double[] thisActions =  new double[] {y1, y2};
					double thisActionsValue = 0;
					for (int j = 0; j < dAndP.length; j++) {
						double[] thisDemands = new double[] { dAndP[j][0], dAndP[j][1] };
						CashStateMulti newState = stateTransition.apply(s, thisActions, thisDemands);
						if (s.getPeriod() < T)
							thisActionsValue += dAndP[j][2] * discountFactor * getExpectedValueV(newState);
						else
							thisActionsValue += dAndP[j][2] * discountFactor * boundFinalCash.apply(newState);
					}
					if (thisActionsValue > bestValue - 0.1) {
						bestValue = thisActionsValue;
						bestAlpha = alpha;
					}
				}
				cacheAlpha.putIfAbsent(initialState, bestAlpha);
			}

			return bestYStar;	
		});
	}
	
	/**
	 * 
	 * @return more detailed optimal decision table
	 */
	public double[][] getOptTableDetail(double[] mean, double[] variance, double[] price, double[] a1, double[] a2){
		Iterator<Map.Entry<CashStateMulti, double[]>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][13]; // revise
		int i = 0;
		Map.Entry<CashStateMulti, double[]> entry;
		while (iterator.hasNext()) {
			entry = iterator.next();
			double[] yStars = cacheYStar.get(entry.getKey());
			double y1 = entry.getValue()[0]; double y2 = entry.getValue()[1];			
			
			double cashConstrained = 0;
			double alpha = 10000;
			double x1 = entry.getKey().iniInventory1; double x2 = entry.getKey().iniInventory2;
			double w = entry.getKey().getIniCash(); 
			double[] yHeads = new double[] {0, 0};
			CashStateMulti s =  entry.getKey();
			if (x1 < yStars[0]+0.1 && x2 < yStars[1]
					&& variCost[0] * (yStars[0] - x1) + variCost[1] * (yStars[1] - x2) < w+0.1) {
				yHeads[0] = yStars[0]; yHeads[1] = yStars[1];
				cashConstrained = 1;
			}
			else if (x1 > yStars[0]-0.1 && x2 > yStars[1]-0.1) {
				yHeads[0] = x1; yHeads[1] = x2;
				cashConstrained = 5;
			}
			else if (x1 > yStars[0]-0.1 && x2 < yStars[1]+0.1) {
				yHeads[0] = x1; yHeads[1] = Math.min(yStars[1], (s.iniCash + x2*variCost[1]) / variCost[1]);
				cashConstrained = 4;
			}
			else if (x1 < yStars[0]+0.1 && x2 > yStars[1]-0.1){
				yHeads[0] = Math.min(yStars[0], (s.iniCash + x1*variCost[0]) / variCost[0]); yHeads[1] = x2;
				cashConstrained = 3;
			}
			else if (x1 < yStars[0]+0.1 && x2 < yStars[1]+0.1
					&& variCost[0] * (yStars[0] - x1) + variCost[1] * (yStars[1] - x2)> s.iniCash-0.1) {
				alpha = cacheAlpha.get(s);
				yHeads[0] = alpha * s.iniCash / variCost[0] + x1; 
				yHeads[1] = (1 - alpha) * s.iniCash / variCost[1] + x2;	
				cashConstrained = 2;
			}				
			int period = s.period;
			double R = s.getIniR(variCost);
			arr[i++] = new double[]{mean[0], mean[1], variance[0], variance[1], s.period, x1, x2, w, price[0], price[1], variCost[0], variCost[1], R, yStars[0], yStars[1], cashConstrained, alpha, yHeads[0], yHeads[1], a1[period-1], a2[period-1]};					
		}
		return arr;
	}

}
