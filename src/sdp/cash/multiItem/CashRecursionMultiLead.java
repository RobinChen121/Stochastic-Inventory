package sdp.cash.multiItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;

import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;

/**
*@author: zhenchen
*@date: Apr 8, 2024, 4:29:02 PM
*@desp: TODO
*
*/

public class CashRecursionMultiLead {
	double discountFactor;	
	int TLength;
	GetPmfMulti Pmf;
	double[][][] pmf;
	Map<CashStateMultiLead, Actions> cacheActions = new ConcurrentSkipListMap<>();	
	Map<CashStateMultiLead, Double> cacheValues = new ConcurrentSkipListMap<>();
	
	Function<CashStateMultiLead, ArrayList<Actions>> buildActionList;
	StateTransitionFunction<CashStateMultiLead, Actions, Demands, CashStateMultiLead> stateTransition;
	ImmediateValueFunction<CashStateMultiLead, Actions, Demands, Double> immediateValue;
	
	public CashRecursionMultiLead(double discountFactor, GetPmfMulti Pmf, Function<CashStateMultiLead, ArrayList<Actions>> buildActionList,
			StateTransitionFunction<CashStateMultiLead, Actions, Demands, CashStateMultiLead> stateTransition, 
			ImmediateValueFunction<CashStateMultiLead, Actions, Demands, Double> immediateValue, int TLength) {
		this.discountFactor = discountFactor;
		this.Pmf = Pmf;
		this.buildActionList = buildActionList;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		this.TLength = TLength;
		
		// sorted map for recorded actions and values 
		Comparator<CashStateMultiLead> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory1() > o2.getIniInventory1() ? 1 : 
				o1.getIniInventory1() == o2.getIniInventory1() ? o1.getIniInventory2() > o2.getIniInventory2() ? 1 :
					o1.getIniInventory2() == o2.getIniInventory2() ? o1.getPreQ1() > o2.getPreQ1() ? 1 :
						o1.getPreQ1() == o2.getPreQ1() ? o1.getPreQ2() > o2.getPreQ2() ? 1 :
							o1.getPreQ2() == o2.getPreQ2() ? o1.iniCash > o2.iniCash  ? 1 :
					o1.iniCash == o2.iniCash ? 0 : -1: -1 : -1 : -1 : -1 : -1;
		this.cacheActions = new ConcurrentSkipListMap<>(keyComparator);
		this.cacheValues = new ConcurrentSkipListMap<>(keyComparator);		
	}
	
	
	public double getExpectedValue(CashStateMultiLead initialState) {
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
				//thisActions = new Actions(9, 8);
				double thisActionsValue = 0;
				for (int j = 0; j < dAndP.length; j++) {
					Demands thisDemands = new Demands((int) dAndP[j][0], (int) dAndP[j][1]);
					thisActionsValue += dAndP[j][2] * immediateValue.apply(s, thisActions, thisDemands);
					if (s.getPeriod()  < TLength) {
						CashStateMultiLead newState = stateTransition.apply(s, thisActions, thisDemands);
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
	
	
	public Actions getAction(CashStateMultiLead state) {
		return cacheActions.get(state);
	}
	
	
	

}


