/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 24, 2019, 11:40:52 AM
 * @Desc: 
 *
 *
 * 
 */
package sdp.cash.multiItem;


import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdistmulti.BiNormalDist;

/**
 * get demands possibilities for two items in a given period
 */
public class GetPmfMulti {
	BiNormalDist[] distributions;
	double truncationQuantile;
	double stepSize;
	
	
	public GetPmfMulti(BiNormalDist[] distributions, double truncationQuantile, double stepSize) {
		this.distributions = distributions;
		this.truncationQuantile = truncationQuantile;
		this.stepSize = stepSize;
	}
	
	public double[][] getPmf(int t){
		
		BiNormalDist tDistributions = distributions[t];

		
		if (t > 2)
			stepSize = 2;
		if (t > 4)
			stepSize = 4;
		
		
		NormalDist distribution1 = new NormalDist(tDistributions.getMu1(), tDistributions.getSigma2());
		NormalDist distribution2 = new NormalDist(tDistributions.getMu2(), tDistributions.getSigma2());
		
		double[] supportLB = new double[2];
		double[] supportUB = new double[2];
		supportLB[0] = (int)distribution1.inverseF(1 - truncationQuantile);
		supportUB[0] = (int)distribution1.inverseF(truncationQuantile);
		supportLB[1] = (int)distribution2.inverseF(1 - truncationQuantile);
		supportUB[1] = (int)distribution2.inverseF(truncationQuantile);
		
		int demandLength1 = (int) ((supportUB[0] - supportLB[0] + 1) / stepSize);
		int demandLength2 = (int) ((supportUB[1] - supportLB[1] + 1) / stepSize);
		
		double[][] pmf = new double[demandLength1 * demandLength2][3];
		int index = 0;
		for (int i = 0; i< demandLength1; i++)
			for (int j = 0; j < demandLength2; j++) {				
				pmf[index][0] = supportLB[0] + i * stepSize;
				pmf[index][1] = supportLB[1] + j * stepSize;
				double probilitySum = truncationQuantile * truncationQuantile;
				pmf[index][2] = (distribution1.cdf(pmf[index][0] + 0.5 * stepSize)
						- distribution1.cdf(pmf[index][0] - 0.5 * stepSize)) * (distribution2.cdf(pmf[index][1] + 0.5 * stepSize)
								- distribution2.cdf(pmf[index][1] - 0.5 * stepSize)) / probilitySum;							
				index++;
			}
		
		return pmf;
	}
}
