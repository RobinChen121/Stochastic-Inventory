package sdp.inventory;

import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.probdist.UniformIntDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 9, 2018---12:48:56 PM
 * @description: get probabilities for different demands in different periods,
 *               and for multi-variate distribution
 */


public class GetPmf {
	Distribution[] distributions;
	Distribution[][] distributionsMulti;
	Distribution distribution;
	double truncationQuantile;
	double stepSize;

	public GetPmf(Distribution[] distributions, double truncationQuantile, double stepSize) {
		this.distributions = distributions;
		this.truncationQuantile = truncationQuantile;
		this.stepSize = stepSize;
	}
	

	public GetPmf(Distribution[][] distributions, double truncationQuantile, double stepSize) {
		this.distributionsMulti = distributions;
		this.truncationQuantile = truncationQuantile;
		this.stepSize = stepSize;
	}
	
	public GetPmf(Distribution distribution, double truncationQuantile, double stepSize) {
		this.truncationQuantile = truncationQuantile;
		this.stepSize = stepSize;
		this.distribution = distribution;
	}
	
	/**
	 * for a single distribution
	 * @return
	 */
	public double[][] getpmfSingleDist(){
		double supportLB;
		double supportUB;
		if (distribution instanceof DiscreteDistributionInt)
			supportLB = 0;
		else {
			supportLB = (int) distribution.inverseF(1 - truncationQuantile);
		}
		supportUB = (int) distribution.inverseF(truncationQuantile);
		int demandLength = (int) ((supportUB - supportLB + 1) / stepSize);
		double[][] pmf = new double[demandLength][];
		
		for (int j = 0; j < demandLength; j++) {
			pmf[j][0] = supportLB + j * stepSize;
		
			if (distribution instanceof DiscreteDistributionInt ||
					distributions[0] instanceof PoissonDist) { // may be something wrong, poisson[] can't be (casted) delivered
				                                                    // but the results are correct
				double probilitySum = distribution.cdf(supportUB) - distribution.cdf(supportLB);
				pmf[j][1] = ((DiscreteDistributionInt) distribution).prob(j) / probilitySum;
			} else {
				double probilitySum = distribution.cdf(supportUB + 0.5 * stepSize)
						- distribution.cdf(supportLB - 0.5 * stepSize);
				pmf[j][1] = (distribution.cdf(pmf[j][0] + 0.5 * stepSize)
						- distribution.cdf(pmf[j][0] - 0.5 * stepSize)) / probilitySum;
			}
		}
		return pmf;
	}
	
	/**
	 * for distributions in a planning horizon
	 * @return
	 */
	public double[][][] getpmf() {
		int T = distributions.length;
		double[] supportLB = new double[T];
		double[] supportUB = new double[T];
		for (int i = 0; i < T; i++) {
			supportLB[i] = (int) distributions[i].inverseF(1 - truncationQuantile);
			if (distributions[0] instanceof DiscreteDistributionInt)
				supportLB[i] = 0;
			supportUB[i] = (int) distributions[i].inverseF(truncationQuantile);
		}

		double[][][] pmf = new double[T][][];
		
		
		
		if (distributions[0] instanceof UniformIntDist) {
			for (int i = 0; i < T; i++) {
				UniformIntDist distribution = (UniformIntDist) distributions[0];
				int demandLength = distribution.getJ() - distribution.getI() + 1;
				pmf[i] = new double[demandLength][];
				int index = 0;
				for (int j = distribution.getXinf(); j <= distribution.getXsup(); j++) {
					pmf[i][index] = new double[2];
					pmf[i][index][0] = j;
					pmf[i][index][1] = distribution.prob(j);
					index++;
				}
			}
			return pmf;
		}
		
		for (int i = 0; i < T; i++) {
			int demandLength = (int) ((supportUB[i] - supportLB[i] + 1) / stepSize);
			pmf[i] = new double[demandLength][];
			// demand values are all integers
			for (int j = 0; j < demandLength; j++) {
				pmf[i][j] = new double[2];
				pmf[i][j][0] = supportLB[i] + j * stepSize;
				if (distributions[0] instanceof DiscreteDistributionInt ||
						distributions[0] instanceof PoissonDist) { // may be something wrong, poisson[] can't be (casted) delivered
					                                                    // but the results are correct
					double probilitySum = distributions[i].cdf(supportUB[i]) - distributions[i].cdf(supportLB[i] - 1);
					pmf[i][j][1] = ((DiscreteDistributionInt) distributions[i]).prob(j) / probilitySum;
				} else {
					double probilitySum = distributions[i].cdf(supportUB[i] + 0.5 * stepSize)
							- distributions[i].cdf(supportLB[i] - 0.5 * stepSize);
					pmf[i][j][1] = (distributions[i].cdf(pmf[i][j][0] + 0.5 * stepSize)
							- distributions[i].cdf(pmf[i][j][0] - 0.5 * stepSize)) / probilitySum;
				}
			}
		}
		return pmf;
	}
	
	
	/**
	* @Description: possibility of demand values for multi-variate distribution
	* @param @return    
	* @return   A 3D matrix
	*/
	public double[][][] getpmfMulti() {
		int T = distributionsMulti.length;
		int N = 2;  //distributionsMulti[0].length; // N is 2 in this code
		
		double[][] supportLB = new double[T][N];
		double[][] supportUB = new double[T][N];
		
		int[][] demandLengths = new int[T][N];
		int[] demandNumt = new int[T];
		for (int t = 0; t < T; t++) {
			int demandLength = 1;
			for (int i =0; i < N; i++) {
				supportLB[t][i] = (int) distributionsMulti[t][i].inverseF(1 - truncationQuantile);
				supportUB[t][i] = (int) distributionsMulti[t][i].inverseF(truncationQuantile);
				demandLengths[t][i]= (int) ((supportUB[t][i] - supportLB[t][i] + 1) / stepSize);
				demandLength *= demandLengths[t][i];
			}
			demandNumt[t] = demandLength;
		}

		double[][][] pmf = new double[T][][];

		for (int t = 0; t < T; t++) {
			pmf[t] = new double[demandNumt[t]][];
			int j = 0;
			int demandSize1 = demandLengths[t][0];
			int demandSize2 = demandLengths[t][1];
			for (int j1 = 0; j1 < demandSize1; j1++) {
				for (int j2 = 0; j2 < demandSize2; j2++) {
					pmf[t][j] = new double[3];
					pmf[t][j][0] = supportLB[t][0] + j1 * stepSize;
					pmf[t][j][1] = supportLB[t][1] + j2 * stepSize;
					int demand1 = (int) pmf[t][j][0];
					int demand2 = (int) pmf[t][j][1];
					double probilitySum1 = distributionsMulti[t][0].cdf(supportUB[t][0]) - distributionsMulti[t][0].cdf(supportLB[t][0]);
					double probilitySum2 = distributionsMulti[t][1].cdf(supportUB[t][1]) - distributionsMulti[t][1].cdf(supportLB[t][1]);
					double lowBound1 = demand1 - 0.5 * stepSize < 0 ? -1 : demand1 - 0.5 * stepSize;
					double prob1 = (distributionsMulti[t][0].cdf(demand1 + 0.5 * stepSize)
							- distributionsMulti[t][0].cdf(lowBound1)) / probilitySum1;
					double lowBound2 = demand2 - 0.5 * stepSize < 0 ? -1 : demand2 - 0.5 * stepSize;
					double prob2 = (distributionsMulti[t][1].cdf(demand2 + 0.5 * stepSize)
							- distributionsMulti[t][1].cdf(lowBound2)) / probilitySum2;				
					pmf[t][j][2] = prob1 * prob2 / (probilitySum1 * probilitySum2);
					j++;
				}
			}	
			
		}
		return pmf;
	}
}
