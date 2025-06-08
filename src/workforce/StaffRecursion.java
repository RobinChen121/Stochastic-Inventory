package workforce;

import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;

import sdp.inventory.State;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.BinomialDist;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;

public class StaffRecursion {
	Map<StaffState, Double> cacheValues = new ConcurrentSkipListMap<>();
	Map<StaffState, Integer> cacheActions = new ConcurrentSkipListMap<>();
	
	Function<StaffState, int[]> getFeasibleAction;
	StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition;
	ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue;
	double[] dimissionRate;
	double[][][][] pmfs;
	int T;
	
	public StaffRecursion(Function<StaffState, int[]> getFeasibleAction, StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition,
			ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue, double[] dimissionRate) {
		this.getFeasibleAction = getFeasibleAction;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		this.dimissionRate = dimissionRate;
		
		Comparator<StaffState> keyComparator = (o1, o2) -> o1.period > o2.period ? 1 : 
			o1.period == o2.period ? o1.iniStaffNum > o2.iniStaffNum ? 1 : 
				o1.iniStaffNum == o2.iniStaffNum ? 0 : -1 : -1;
		this.cacheActions = new ConcurrentSkipListMap<>(keyComparator);
		this.cacheValues = new ConcurrentSkipListMap<>(keyComparator);
	}
	
	public StaffRecursion(Function<StaffState, int[]> getFeasibleAction, StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition,
			ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue, double[][][][] pmf, int T) {
		this.getFeasibleAction = getFeasibleAction;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		this.pmfs = pmf;
		this.T = T;
		
		Comparator<StaffState> keyComparator = (o1, o2) -> o1.period > o2.period ? 1 : 
			o1.period == o2.period ? o1.iniStaffNum > o2.iniStaffNum ? 1 : 
				o1.iniStaffNum == o2.iniStaffNum ? 0 : -1 : -1;
		this.cacheActions = new ConcurrentSkipListMap<>(keyComparator);
		this.cacheValues = new ConcurrentSkipListMap<>(keyComparator);
	}
	
	public StateTransitionFunction<StaffState, Integer, Integer, StaffState> getStateTransitionFunction(){
		return stateTransition;
	}
	
	public ImmediateValueFunction<StaffState, Integer, Integer, Double> getImmediateValueFunction(){
		return immediateValue;
	}
	
	public double[][] getPmf(Distribution distribution, int n){
		int demandLength = n + 1;
		double[][] pmf = new double[demandLength][2];
		
		for (int j = 0; j < demandLength; j++) {
			pmf[j][0] = j;		
			pmf[j][1] = ((DiscreteDistributionInt) distribution).prob(j);
		}
		return pmf;
	}
	
	/**
	 * binomial demand is related with y
	 * @date ：2022 Aug 22 16:21:10
	 * @param state
	 * @return
	 */
	public double getExpectedValue(StaffState state) {
		return this.cacheValues.computeIfAbsent(state, s -> {	
			int[] feasibleActions = getFeasibleAction.apply(state);
			int iniStaffNum = state.iniStaffNum;
			int t = state.period - 1;
			
			int bestHireQty = 0;
			double[] QValues = new double[feasibleActions.length];
			double val = Double.MAX_VALUE;
			for (int i = 0; i < feasibleActions.length; i++) {
				int orderQty = feasibleActions[i];
				int hireUpTo = iniStaffNum + orderQty;
				if (hireUpTo >= pmfs[t].length - 1) // important iniStaffNum + orderQty >= pmfs[t].length
					hireUpTo = pmfs[t].length - 1;
				double[][] pmf = pmfs[t][hireUpTo];
				double thisQValue = 0;
				for (int j = 0; j < pmf.length; j++) {
					int demand = (int) pmf[j][0];
					double thisValue = immediateValue.apply(s, orderQty, demand);
					thisQValue += pmf[j][1] * thisValue;

					if (s.period < T) {
						StaffState newState = stateTransition.apply(s, orderQty, demand);
						thisQValue += pmf[j][1] * getExpectedValue(newState);
					}
				}
				if (t == 0 && i >= 119) {
					int a = 0;
				}
			        
				QValues[i] = thisQValue;
				if (QValues[i] < val) {
					val = QValues[i];
					bestHireQty = orderQty;
				}
			}
			
			this.cacheActions.putIfAbsent(s, bestHireQty);
			return val;
		});
	}
		
	/**
	 * for drawing G(y), demand is related with y
	 * @date ：2022 Aug 22 17:00:00
	 * @param state
	 * @param hireUpStaffNum
	 * @return
	 */
	public double getExpectedValue(StaffState state, int hireUpStaffNum) {
		return this.cacheValues.computeIfAbsent(state, s -> {	
			int[] feasibleActions = getFeasibleAction.apply(s);
			int t = state.period - 1;
					
			int bestHireQty = 0;
			double[] QValues = new double[feasibleActions.length];
			double val = Double.MAX_VALUE;
			for (int i = 0; i < feasibleActions.length; i++) {
				int orderQty = feasibleActions[i];
				double[][] pmf;
				int hireUpTo = s.iniStaffNum + orderQty;
				if (s.iniStaffNum + orderQty >= pmfs[t].length - 1) 
					hireUpTo = pmfs[t].length - 1;
				pmf = pmfs[t][hireUpTo];
//				if (s.period == 1) {
//					if (s.iniStaffNum  >= pmfs[t].length)
//						continue;
//					pmf = pmfs[t][s.iniStaffNum];
//				}
//				else {
//					if (s.iniStaffNum + orderQty >= pmfs[t].length) 
//						continue;
//					pmf = pmfs[t][s.iniStaffNum + orderQty];
//				}
				
				double thisQValue = 0;								
				for (int j = 0; j < pmf.length; j++) {
					int demand = (int) pmf[j][0];
					double thisValue = immediateValue.apply(s, orderQty, demand);
					thisQValue += pmf[j][1] * thisValue;
					
					if (s.period < T) {
						StaffState newState = stateTransition.apply(s, orderQty, demand);
						thisQValue += pmf[j][1] * getExpectedValue(newState);
					}
				}
				QValues[i] = thisQValue;
				if (QValues[i] < val) {
					val = QValues[i];
					bestHireQty = orderQty;
				}
			}			
			this.cacheActions.putIfAbsent(s, bestHireQty);
			return val;
		});
	}
	
	/**
	 * for drawing G(y), demand is related with x, x is realStaffNum
	 * @date ：2022 Aug 22 17:00:00
	 * @param state
	 * @param hireUpStaffNum
	 * @return
	 */
	public double getExpectedValue2(StaffState state, int realStaffNum) {
		return this.cacheValues.computeIfAbsent(state, s -> {	
			int[] feasibleActions = getFeasibleAction.apply(s);
			int iniStaffNum = s.period == 1 ? realStaffNum : s.iniStaffNum;
			int t = state.period - 1;
			double[][] pmf = {{0, 1}};
			if (s.period == 1) {
				if (realStaffNum > 0) {
					Distribution distribution = new BinomialDist(realStaffNum , dimissionRate[t]);
					pmf = getPmf(distribution, realStaffNum);
				}
			}
			else {
				if (iniStaffNum > 0) {
				Distribution distribution = new BinomialDist(iniStaffNum, dimissionRate[t]);
				pmf = getPmf(distribution, iniStaffNum);
				}
			}		
			
			int bestHireQty = 0;
			double[] QValues = new double[feasibleActions.length];
			double val = Double.MAX_VALUE;
			for (int i = 0; i < feasibleActions.length; i++) {
				int orderQty = feasibleActions[i];
				double thisQValue = 0;								
				for (int j = 0; j < pmf.length; j++) {
					double thisValue = immediateValue.apply(s, orderQty, (int)pmf[j][0]);
					thisQValue += pmf[j][1] * thisValue;

					if (s.period < dimissionRate.length) {
						StaffState newState = stateTransition.apply(s, orderQty, (int) pmf[j][0]);
						thisQValue += pmf[j][1] * getExpectedValue(newState);
					}
				}
				QValues[i] = thisQValue;
				if (QValues[i] < val) {
					val = QValues[i];
					bestHireQty = orderQty;
				}
			}
			
			this.cacheActions.putIfAbsent(s, bestHireQty);
			return val;
		});
	}
	
	/**
	 * expected value when not hiring in the first period
	 * @param state
	 * @return
	 */
	public double getExpectedValueNoHireFirst(StaffState state) {
		return this.cacheValues.computeIfAbsent(state, s -> {	
			int[] feasibleActions = getFeasibleAction.apply(state);
			int iniStaffNum = state.iniStaffNum;
			int t = state.period - 1;
			Distribution distribution = new BinomialDist(iniStaffNum, dimissionRate[t]);
			double[][] pmf = getPmf(distribution, iniStaffNum);
			
			int bestHireQty = 0;
			double[] QValues = new double[feasibleActions.length];
			double val = Double.MAX_VALUE;
			for (int i = 0; i < feasibleActions.length; i++) {
				int orderQty = t == 0 ? 0 : feasibleActions[i];
				double thisQValue = 0;								
				for (int j = 0; j < pmf.length; j++) {
					double thisValue = immediateValue.apply(s, orderQty, (int)pmf[j][0]);
					thisQValue += pmf[j][1] * thisValue;
					if (s.period < dimissionRate.length) {
						StaffState newState = stateTransition.apply(s, orderQty, (int) pmf[j][0]);
						thisQValue += pmf[j][1] * getExpectedValue(newState);
					}
				}
				QValues[i] = thisQValue;
				if (QValues[i] < val) {
					val = QValues[i];
					bestHireQty = orderQty;
				}
			}

			return val;
		});
	}
	
	public int getAction(StaffState state) {
		return cacheActions.get(state);
	}
	
	/**
	 * 
	 * @return optimal decision table of SDP
	 */
	public double[][] getOptTable(){
		Iterator<Map.Entry<StaffState, Integer>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][3];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<StaffState, Integer> entry = iterator.next();
			arr[i++] = new double[]{entry.getKey().period, entry.getKey().iniStaffNum, entry.getValue()};
		}
		return arr;
	}
	
	/**
	 * binomial demand is related with x
	 * @date ：2022 Aug 22 16:21:10
	 * @param state
	 * @return
	 */
	public double getExpectedValue2(StaffState state) {
		return this.cacheValues.computeIfAbsent(state, s -> {	
			int[] feasibleActions = getFeasibleAction.apply(state);
			int iniStaffNum = state.iniStaffNum;
			int t = state.period - 1;
			double[][] pmf = {{0, 1}};
			if (iniStaffNum > 0) {
				Distribution distribution = new BinomialDist(iniStaffNum, dimissionRate[t]);
				pmf = getPmf(distribution, iniStaffNum);
			}
			
			int bestHireQty = 0;
			double[] QValues = new double[feasibleActions.length];
			double val = Double.MAX_VALUE;
			for (int i = 0; i < feasibleActions.length; i++) {
				int orderQty = feasibleActions[i];
				double thisQValue = 0;								
				for (int j = 0; j < pmf.length; j++) {
					double thisValue = immediateValue.apply(s, orderQty, (int)pmf[j][0]);
					try {
						thisQValue += pmf[j][1] * thisValue;
					} catch (Exception e) {
						System.out.println();
					}

					if (s.period < dimissionRate.length) {
						StaffState newState = stateTransition.apply(s, orderQty, (int) pmf[j][0]);
						thisQValue += pmf[j][1] * getExpectedValue2(newState);
					}
				}
				QValues[i] = thisQValue;
				if (QValues[i] < val) {
					val = QValues[i];
					bestHireQty = orderQty;
				}
			}
			
			this.cacheActions.putIfAbsent(s, bestHireQty);
			return val;
		});
	}

}
