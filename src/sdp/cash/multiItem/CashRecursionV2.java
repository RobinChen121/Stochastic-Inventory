package sdp.cash.multiItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;

import sdp.inventory.FinalCash.BoundaryFuncton;
import sdp.inventory.StateTransition.StateTransitionFunction;

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
	Map<CashStateMulti, double[]> cacheActions = new TreeMap<>();	
	Map<CashStateMulti, double[]> cacheYStar = new TreeMap<>();	
	Map<CashStateMulti, Double> cacheAlpha = new TreeMap<>();	
	Map<CashStateMulti, Double> cacheValuesV = new TreeMap<>();
	Map<CashStateMulti, Double> cacheValuesPai = new TreeMap<>();
	
	Function<CashStateMulti, ArrayList<double[]>> buildActionListV;
	Function<CashStateMulti, ArrayList<double[]>> buildActionListPai;
	StateTransitionFunction<CashStateMulti, double[], double[], CashStateMulti> stateTransition;
	BoundaryFuncton<CashStateMulti, Double> boundFinalCash;
	
	public CashRecursionV2(double discountFactor, GetPmfMulti Pmf, Function<CashStateMulti, ArrayList<double[]>> buildActionListV,
			Function<CashStateMulti, ArrayList<double[]>> buildActionListPai, StateTransitionFunction<CashStateMulti, double[], double[],CashStateMulti> stateTransition, 
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
					
		this.cacheActions = new TreeMap<>(keyComparator1);
		this.cacheValuesV = new TreeMap<>(keyComparator1);		
		this.cacheValuesPai = new TreeMap<>(keyComparator1);	
		this.cacheYStar = new TreeMap<>(keyComparator1);
		this.cacheAlpha = new TreeMap<>(keyComparator1);		
	}
	
	
	public double getExpectedValuePai(CashStateMulti initialState, double[] actions) {
		return this.cacheValuesPai.computeIfAbsent(initialState, s -> {
			int n = s.getPeriod();
			double[][] dAndP = Pmf.getPmf(n - 1); // demandAndPossibility
			double expectValue = 0;
			for (int j = 0; j < dAndP.length; j++) {
				double[] thisDemands = new double[] {dAndP[j][0],  dAndP[j][1]};
				CashStateMulti newState = stateTransition.apply(s, actions, thisDemands);
				//double thisProfit = immediateValue.apply(newState, thisDemands);
				double thisDemandValue = getExpectedValueV(newState);
				expectValue += dAndP[j][2] * thisDemandValue;
			}	
			return expectValue;
		});
	}

	
	public double getExpectedValueV(CashStateMulti initialState) {
		return this.cacheValuesV.computeIfAbsent(initialState, s -> {
			ArrayList<double[]> yHeads = buildActionListV.apply(s);
			double val = -Double.MAX_VALUE;
			double[] bestYs = new double[] {initialState.getIniInventory1(), initialState.getIniInventory2()};
			if (initialState.getPeriod() <= T) {
				for (int i = 0; i < yHeads.size(); i++) {
					double[] thisActions = yHeads.get(i);
					CashStateMulti thisState = new CashStateMulti(s.getPeriod(), s.iniInventory1, s.iniInventory2, s.iniCash);
					double thisActionsValue = getExpectedValuePai(thisState, thisActions);

					if (thisActionsValue > val + 0.1) {
						val = thisActionsValue;
						bestYs = thisActions;
					}
				}
				getYStar(initialState); // for cacheing y stars
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
	* @return optimal y1, y2  for state (x1, x2, w)
	*/
	public double[] getAction(CashStateMulti state) {
		return cacheActions.get(state);
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
		double[] bestYs = new double[] {0, 0};
		for (int i = 0; i < actions.size(); i++) {
			double[] thisActions = actions.get(i);
			double thisActionsValue = getExpectedValuePai(initialState, thisActions);
			
			if (thisActionsValue > val + 0.01) {
				val = thisActionsValue;
				bestYs = thisActions;
			}					
		}

		if (variCost[0] * (bestYs[0] - s.iniInventory1) + variCost[1] * (bestYs[1] - s.iniInventory2) >= s.iniCash + 0.1) {
			getAlpha(s); // revise to save computation time
		}
			
		return bestYs;			
		});
	}
	

	public double getAlpha(CashStateMulti initialState) {
		return this.cacheAlpha.computeIfAbsent(initialState, s -> {
			double bestAlpha = 0;
			double bestValue = -Double.MAX_VALUE;
			for (double alpha = 0; alpha <= 1; alpha = alpha + 0.01) {  // stepsize of alpha
				double y1 = alpha * s.iniCash / variCost[0] + s.iniInventory1;
				double y2 = (1 - alpha) * s.iniCash / variCost[1] + s.iniInventory2;
				double[] actions =  new double[] {y1, y2};
				double expectValue = getExpectedValuePai(initialState, actions);
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
				yHeads[0] = Math.min(yStars[0], (s.iniCash -  x1*variCost[0]) / variCost[0]); yHeads[1] = x2;
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
