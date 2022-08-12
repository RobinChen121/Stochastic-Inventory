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
	double truncQuantile;
	
	public StaffRecursion(Function<StaffState, int[]> getFeasibleAction, StateTransitionFunction<StaffState, Integer, Integer, StaffState> stateTransition,
			ImmediateValueFunction<StaffState, Integer, Integer, Double> immediateValue, double[] dimissionRate, double truncQuantile) {
		this.getFeasibleAction = getFeasibleAction;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		this.dimissionRate = dimissionRate;
		this.truncQuantile = truncQuantile;
		
		Comparator<StaffState> keyComparator = (o1, o2) -> o1.period > o2.period ? 1 : 
			o1.period == o2.period ? o1.iniStaffNum > o2.iniStaffNum ? 1 : 
				o1.iniStaffNum == o2.iniStaffNum ? 0 : -1 : -1;
		this.cacheActions = new ConcurrentSkipListMap<>(keyComparator);
		this.cacheValues = new ConcurrentSkipListMap<>(keyComparator);
	}
	
	public double[][] getPmf(Distribution distribution){
		double supportLB;
		double supportUB;
		
		supportLB = 0;
		supportUB = (int) distribution.inverseF(truncQuantile);
		int demandLength = (int) (supportUB - supportLB + 1);
		double[][] pmf = new double[demandLength][2];
		
		for (int j = 0; j < demandLength; j++) {
			pmf[j][0] = supportLB + j;		
			double probilitySum = distribution.cdf(supportUB) - distribution.cdf(supportLB);
			pmf[j][1] = ((DiscreteDistributionInt) distribution).prob(j) / probilitySum;
		}
		return pmf;
	}
	
	public double getExpectedValue(StaffState state) {
		return this.cacheValues.computeIfAbsent(state, s -> {	
			int[] feasibleActions = getFeasibleAction.apply(state);
			int iniStaffNum = state.iniStaffNum;
			int t = state.period - 1;
			double[][] pmf = {{0, 1}};; 
			if (iniStaffNum > 0) {
				Distribution distribution = new BinomialDist(iniStaffNum, dimissionRate[t]);
				pmf = getPmf(distribution);
			}
			else {
				
			}
			
			int bestHireQty = 0;
			double[] QValues = new double[feasibleActions.length];
			double val = Double.MAX_VALUE;
//			if (t == 1 && Math.abs(state.iniStaffNum - 93) < 0.1)
//				System.out.println(state.iniStaffNum);
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
			double[][] pmf = getPmf(distribution);
			
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

}
