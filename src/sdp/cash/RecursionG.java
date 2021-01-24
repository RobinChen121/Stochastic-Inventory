
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
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Apr 17, 2020, 8:43:09 PM
 * @Desc: testing the algorithm proposed by Chao et al. (2008)
 *
 *
 * 
 */

public class RecursionG {
	
	Map<StateP, Double> cachePeriodBestY = new TreeMap<>();	
	Map<StateY, Double> cacheGValues = new TreeMap<>();
	
	
	double[][][] pmf;
	double tOptY[]; // optimal Y in each period
	
	double price;
	double variCost;
	double depositeRate;
	double salvageValue;
	
	
	Distribution[] distributions;
	
	double aNStar;
	
	public RecursionG(double[][][] pmf, Distribution[] distributions,
			double price, double variCost, double depositeRate, double salvageValue) {
		this.pmf = pmf;
		this.aNStar = aNStar;
		Comparator<StateY> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniY() > o2.getIniY() ? 1 : 
				o1.getIniY() == o2.getIniY() ? 0 : -1 : -1;
 		Comparator<StateP> keyComparator2 = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? 0 : -1;
		this.cachePeriodBestY = new TreeMap<>(keyComparator2);
		this.cacheGValues = new TreeMap<>(keyComparator);
		this.tOptY = new double[pmf.length];
		Arrays.fill(tOptY, -100); // initialize tOptY
		this.distributions = distributions;
		this.price = price;
		this.variCost = variCost;
		this.depositeRate = depositeRate;
		this.salvageValue = salvageValue;
		// avoid infinity
		if (salvageValue < variCost)
			this.aNStar = distributions[distributions.length - 1].inverseF((price - (1 + depositeRate) * variCost) / (price - salvageValue));
		else {
			this.aNStar = distributions[distributions.length - 1].inverseF(0.999);
		}
		
	}
	

	/**
	 * 
	 * @return optimal y* in each period
	 * @date: Apr 23, 2020
	 */	
	public double[] getOptY() {
		int T = distributions.length;
		tOptY[T - 1] = aNStar;
		for (int t = T - 2; t >= 0; t--) {
			StateP newStateP = new StateP(t + 1);
			tOptY[t] = getAStar(newStateP);
		}
		return tOptY;
	}
 
	/**
	 * the two functions call each other: getAStar and G
	 * @param state
	 * @return
	 * @date: Apr 23, 2020
	 */
	public double G(StateY state) {
		return this.cacheGValues.computeIfAbsent(state, s -> {
			int T = distributions.length;
			int n = s.getPeriod();
			double y = s.getIniY();
			double[][] dAndP = pmf[n - 1]; 
			double expectValue = 0;
			if (n == T) {
				for (int j = 0; j < dAndP.length; j++) {
					double thisDValue = (price - variCost) * Math.min(dAndP[j][0], y) - depositeRate * variCost * y
							+ (salvageValue - variCost) * Math.max(y - dAndP[j][0], 0);
					expectValue += dAndP[j][1] * thisDValue;
				}
				return expectValue;
			}
			else {
				StateP newStateP = new StateP(s.getPeriod() + 1);
				double nextAStar = getAStar(newStateP);
				for (int j = 0; j < dAndP.length; j++) {		
					StateY newStateY = new StateY(n + 1, Math.max(nextAStar, Math.max(y - dAndP[j][0], 0)));
					double thisDValue = Math.pow(1 + depositeRate, T - n) * ((price - variCost) * Math.min(dAndP[j][0], y) - depositeRate * variCost * y) 
							+ G(newStateY);
					expectValue += dAndP[j][1] * thisDValue;
				}
				return expectValue;			
			}
		});
	}
	
	/**
	 * the two functions call each other: getAStar and G
	 * @param state
	 * @return a* for a period
	 */
	public double getAStar(StateP state) {
		return this.cachePeriodBestY.computeIfAbsent(state, s ->{
			int T = distributions.length;
			int n = s.getPeriod();
			if (n == T)
				return aNStar;
			int maxY = 200;
			double optY = 0;
			double optYValue = -1000;
			for (double y  = 0; y < maxY; y = y + 0.1) {  // step size
				StateY newStateY = new StateY(s.getPeriod(), y);
				double thisYValue = G(newStateY); // optimal y in the next period
//				if (s.period == 1 && (int) y == 21 || (int) y == 14)
//					System.out.println(thisYValue);
//				if (s.period == 1 && (int) y == 21 || (int) y == 9)
//					System.out.println(thisYValue);
				if (thisYValue - optYValue > 0.01) {
					optYValue = thisYValue;
					optY = y;
				}
			}
			return optY;
		});	
	}
	
	
	/**
	 * @param a* in each period
	 * @return simulate results for the policy proposed by Chao (2008)
	 */
	public double simulateAStar(double[] optY, int sampleNum) {
		Sampling.resetStartStream();
		double[][] samples = Sampling.generateLHSamples(distributions, sampleNum);
		
		double[] simValues = new double[samples.length];
		for (int i = 0; i < samples.length; i++) {}
		return 0;
	}
}
