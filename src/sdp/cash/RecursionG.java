/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Apr 17, 2020, 8:43:09 PM
 * @Desc: testing the algorithm proposed by Chao et al. (2008)
 *
 *
 * 
 */
package sdp.cash;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.BiFunction;
import java.util.function.Function;

import sdp.inventory.State;
import sdp.inventory.StateTransition;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;

public class RecursionG {
	
	Map<StateP, Double> cachePeriodBestY = new TreeMap<>();	
	Map<StateY, Double> cacheValues = new TreeMap<>();
	
	
	double[][][] pmf;
	double tOptY[]; // optimal Y in each period
	
	Function<StateY, double[]> getFeasibleActions;
	StateTransitionFunction<StateY, Double, Double, StateY> stateTransition;
	BiFunction<StateY, Double, Double> immediateValue;
	
	Distribution[] distributions;
	
	double aNStar;
	
	public RecursionG(double[][][] pmf, double aNStar,
			Function<StateY, double[]> getFeasibleAction,
			StateTransitionFunction<StateY, Double, Double, StateY> stateTransition,
			BiFunction<StateY, Double, Double> immediateValue) {
		this.pmf = pmf;
		this.aNStar = aNStar;
		this.getFeasibleActions = getFeasibleAction;
		this.stateTransition = stateTransition;
		this.immediateValue = immediateValue;
		Comparator<StateY> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniY() > o2.getIniY() ? 1 : 
				o1.getIniY() == o2.getIniY() ? 0 : -1 : -1;
		this.cachePeriodBestY = new TreeMap<>();
		this.cacheValues = new TreeMap<>(keyComparator);
		this.tOptY = new double[pmf.length + 1];
		Arrays.fill(tOptY, -100); // initialize tOptY
	}
	
	public StateTransitionFunction<StateY, Double, Double, StateY> getStateTransitionFunction(){
		return stateTransition;
	}
	
	public double getExpectedValue(StateY state) {
		return this.cacheValues.computeIfAbsent(state, s -> {
			this.cachePeriodBestY.putIfAbsent(new StateP(pmf.length), aNStar);
			double[] feasibleActions = getFeasibleActions.apply(state);
			double[][] dAndP = pmf[s.getPeriod() - 1]; // demandAndPossibility
			double[] YValues = new double[feasibleActions.length];			
			double bestY = 0;
			double val = -Double.MAX_VALUE;
			
			double thisYValue = 0;
			double nextOptY = 0;
			
			
			for (int i  = 0; i < feasibleActions.length; i++) {
				double yAction = feasibleActions[i]; 
				double thisActionValue = 0;
				StateY yStateY = new StateY(s.getPeriod(), yAction);
				for (int j = 0; j < dAndP.length; j++) {					
					thisActionValue += dAndP[j][1] * immediateValue.apply(yStateY, dAndP[j][0]);
					if (yStateY.getPeriod() < pmf.length) {	
						nextOptY = getNextOptY(new StateP(yStateY.getPeriod()));
						StateY newState = stateTransition.apply(yStateY, nextOptY, dAndP[j][0]);
						thisActionValue += dAndP[j][1] * getExpectedValue(newState);
					}		
				}
				if (Math.abs(yStateY.getIniY() - yAction) < 0.1)
					thisYValue = thisActionValue;
				YValues[i] = thisActionValue;
				if (YValues[i] > val) {
					val = YValues[i];
					bestY = yAction;
				}				
			}		
			this.cachePeriodBestY.putIfAbsent(new StateP(s.getPeriod()), bestY);
			return thisYValue;
		});
	}
	
	
	public double getNextOptY(StateP state) {
		double nextOptY = 0;
		if (state.getPeriod() == pmf.length - 1)
			nextOptY = aNStar;
		else {
			nextOptY = cachePeriodBestY.get(new StateP(state.getPeriod() + 1));
		}
		
		return nextOptY;
	}

	public double[] getOptY() {
		Iterator<Map.Entry<StateP, Double>> iterator = this.cachePeriodBestY.entrySet().iterator();
		while (iterator.hasNext()) {
			Map.Entry<StateP, Double> entry = iterator.next();
			tOptY[entry.getKey().getPeriod() - 1] = entry.getValue();
		}
		return tOptY;
	}
 
	
	
	
	

}
