package sdp.inventory;

import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.Distribution;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 9, 2018---12:48:56 PM
 * @description: get probabilities for different demands in different periods
 */

public class GetPmf {
	Distribution[] distributions;
	double truncationQuantile;
	double stepSize;

	public GetPmf(Distribution[] distributions, double truncationQuantile, double stepSize) {
		this.distributions = distributions;
		this.truncationQuantile = truncationQuantile;
		this.stepSize = stepSize;
	}

	public double[][][] getpmf() {
		int T = distributions.length;
		double[] supportLB = new double[T];
		double[] supportUB = new double[T];
		for (int i = 0; i < T; i++) {
			supportLB[i] = (int) distributions[i].inverseF(1 - truncationQuantile);
			supportUB[i] = (int) distributions[i].inverseF(truncationQuantile);
		}

		double[][][] pmf = new double[T][][];

		for (int i = 0; i < T; i++) {
			int demandLength = (int) ((supportUB[i] - supportLB[i] + 1) / stepSize);
			pmf[i] = new double[demandLength][];
			// demand values are all integers
			for (int j = 0; j < demandLength; j++) {
				pmf[i][j] = new double[2];
				pmf[i][j][0] = supportLB[i] + j * stepSize;
				if (distributions[0] instanceof DiscreteDistribution) {
					double probilitySum = distributions[i].cdf(supportUB[i]) - distributions[i].cdf(supportLB[i] - 1);
					pmf[i][j][1] = ((DiscreteDistribution) distributions[i]).prob(j) / probilitySum;
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
}
