package sdp.cash.multiItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;

import sdp.inventory.FinalCash.BoundaryFuncton;
import sdp.inventory.StateTransition.StateTransitionFunction;

/**
 * @author chen
 * @email: 15011074486@163.com
 * @Date: 2021 Feb 27, 23:19:10
 * @Description: Test merging recursion function for V and Pai 
 * 
 */
public class CashRecursionVTest {
	
	double discountFactor;	
	int T;
	GetPmfMulti Pmf;
	double[][][] pmf;
	double[] variCost; 
	Map<CashStateMulti, double[]> cacheActions = new TreeMap<>();	
	Map<CashStateMulti, Double> cacheValuesV = new TreeMap<>();
	
	Function<CashStateMulti, ArrayList<double[]>> buildActionListV;
	StateTransitionFunction<CashStateMulti, double[], double[], CashStateMulti> stateTransition;
	BoundaryFuncton<CashStateMulti, Double> boundFinalCash;
	
	public CashRecursionVTest(double discountFactor, GetPmfMulti Pmf, Function<CashStateMulti, ArrayList<double[]>> buildActionListV,
			 StateTransitionFunction<CashStateMulti, double[], double[], CashStateMulti> stateTransition, 
			BoundaryFuncton<CashStateMulti, Double> boundFinalCash, int T, double[] variCost) {
		this.discountFactor = discountFactor;
		this.Pmf = Pmf;
		this.buildActionListV = buildActionListV;
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
	}
	
	public double getExpectedValueV(CashStateMulti initialState) {
		return this.cacheValuesV.computeIfAbsent(initialState, s -> {
			double val = -Double.MAX_VALUE;
			ArrayList<double[]> yHeads = buildActionListV.apply(s);
			double[][] dAndP = Pmf.getPmf(s.getPeriod() - 1);
			double[] bestYs = new double[] { s.iniInventory1, s.iniInventory2 };

			for (int i = 0; i < yHeads.size(); i++) {
				double[] thisActions = yHeads.get(i);
				double thisActionsValue = 0;
				for (int j = 0; j < dAndP.length; j++) {
					double[] thisDemands = new double[] { dAndP[j][0], dAndP[j][1] };
					CashStateMulti newState = stateTransition.apply(s, thisActions, thisDemands);
					if (s.getPeriod() < T)
						thisActionsValue += dAndP[j][2] * discountFactor * getExpectedValueV(newState);
					else
						thisActionsValue += dAndP[j][2] * discountFactor * boundFinalCash.apply(newState);
				}
				if (thisActionsValue > val) {
					val = thisActionsValue;
					bestYs = thisActions;
				}
			}
			this.cacheActions.putIfAbsent(s, bestYs);	
			return val;
		});
	}
	
	public double[] getAction(CashStateMulti state) {
		return cacheActions.get(state);
	}
	

}
