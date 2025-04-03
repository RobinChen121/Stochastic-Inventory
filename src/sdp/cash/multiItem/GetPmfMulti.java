/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 24, 2019, 11:40:52 AM
 * @Desc: for binormal distribution only
 *
 *
 * 
 */
package sdp.cash.multiItem;



import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.probdist.UniformIntDist;
import umontreal.ssj.probdistmulti.BiNormalDist;

/**
 * get demands possibilities for two items in a given period
 */
public class GetPmfMulti {
	Distribution[][] distributionGeneral;
	double truncationQuantile;
	double stepSize;	
	
	
	
	public GetPmfMulti(Distribution[][] distributions, double truncationQuantile, double stepSize) {
		this.distributionGeneral = distributions;
		this.truncationQuantile = truncationQuantile;
		this.stepSize = stepSize;
	}
	
	
	
	public double[][] getPmf(int t){
//		stepSize = 1;
//		if (t > 4)
//			stepSize = 4;
		
		if(distributionGeneral[0][t] instanceof NormalDist) {
			NormalDist distribution1 = new NormalDist(distributionGeneral[0][t].getMean(), distributionGeneral[0][t].getStandardDeviation());
			NormalDist distribution2 = new NormalDist(distributionGeneral[1][t].getMean(), distributionGeneral[1][t].getStandardDeviation());
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
			double probilitySum = (2 * truncationQuantile - 1) * (2 * truncationQuantile - 1);
			for (int i = 0; i< demandLength1; i++)
				for (int j = 0; j < demandLength2; j++) {				
					pmf[index][0] = supportLB[0] + i * stepSize;
					pmf[index][1] = supportLB[1] + j * stepSize;
					pmf[index][2] = (distribution1.cdf(pmf[index][0] + 0.5 * stepSize)
							- distribution1.cdf(pmf[index][0] - 0.5 * stepSize)) * (distribution2.cdf(pmf[index][1] + 0.5 * stepSize)
									- distribution2.cdf(pmf[index][1] - 0.5 * stepSize)) / probilitySum;							
					index++;
				}	
			return pmf;
		}
		
		if(distributionGeneral[0][t] instanceof GammaDist) {
			
			double mean1 = distributionGeneral[0][t].getMean(); double mean2 = distributionGeneral[1][t].getMean(); 
			double variance1 = distributionGeneral[0][t].getVariance(); double variance2 = distributionGeneral[1][t].getVariance();
			
			double scale1 =  mean1 / variance1;
			double shape1 = mean1 * scale1;
			double scale2 = mean2 / variance2;
			double shape2 = mean2 * scale2;
			
			GammaDist distribution1 = new GammaDist(shape1, scale1);
			GammaDist distribution2 = new GammaDist(shape2, scale2);
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
			double probilitySum = (2 * truncationQuantile - 1) * (2 * truncationQuantile - 1);
			for (int i = 0; i< demandLength1; i++)
				for (int j = 0; j < demandLength2; j++) {				
					pmf[index][0] = supportLB[0] + i * stepSize;
					pmf[index][1] = supportLB[1] + j * stepSize;
					pmf[index][2] = (distribution1.cdf(pmf[index][0] + 0.5 * stepSize)
							- distribution1.cdf(pmf[index][0] - 0.5 * stepSize)) * (distribution2.cdf(pmf[index][1] + 0.5 * stepSize)
									- distribution2.cdf(pmf[index][1] - 0.5 * stepSize)) / probilitySum;							
					index++;
				}	
			return pmf;
		}
		
		if(distributionGeneral[0][t] instanceof PoissonDist) {
			PoissonDist distribution1 = new PoissonDist(distributionGeneral[0][t].getMean());
			PoissonDist distribution2 = new PoissonDist(distributionGeneral[1][t].getMean());
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
			double probilitySum = (2 * truncationQuantile - 1) * (2 * truncationQuantile - 1);
			for (int i = 0; i< demandLength1; i++)
				for (int j = 0; j < demandLength2; j++) {				
					pmf[index][0] = supportLB[0] + i * stepSize;
					pmf[index][1] = supportLB[1] + j * stepSize;
					int demand1 = (int)pmf[index][0];
					int demand2 = (int)pmf[index][1];
					pmf[index][2] = distribution1.prob(demand1) * distribution2.prob(demand2) / probilitySum;
					index++;
				}	
			
//			double pSum = 0;
//			for (int i = 0; i < pmf.length; i++) {
//				pSum += pmf[i][2];
//			}
//			System.out.println(pSum);
			
			return pmf;
		}
		
		if(distributionGeneral[0][t] instanceof UniformIntDist) {
			UniformIntDist distribution1 = (UniformIntDist) distributionGeneral[0][t];
			UniformIntDist distribution2 = (UniformIntDist) distributionGeneral[1][t];
			int demandLength1 = distribution1.getJ() - distribution1.getI() + 1;
			int demandLength2 = distribution2.getJ() - distribution2.getI() + 1;
			double[][] pmf = new double[demandLength1 * demandLength2][3];
			int index = 0;
			for (int i = distribution1.getXinf(); i<= distribution1.getXsup(); i++)
				for (int j = distribution2.getXinf(); j <= distribution2.getXsup(); j++) {	
					pmf[index][0] = i;
					pmf[index][1] = j;
					pmf[index][2] = distribution1.prob(i) * distribution2.prob(j);
					index++;
				}
			return pmf;
		}
	
		if(distributionGeneral[0][t] instanceof DiscreteDistribution) {
			DiscreteDistribution distribution1 = (DiscreteDistribution) distributionGeneral[0][t];
			DiscreteDistribution distribution2 = (DiscreteDistribution) distributionGeneral[1][t];
			int demandLength1 = distribution1.getN();
			int demandLength2 = distribution2.getN();
			double[][] pmf = new double[demandLength1 * demandLength2][3];
			int index = 0;
			for (int i = 0; i < demandLength1; i++)
				for (int j = 0; j < demandLength2; j++) {	
					pmf[index][0] = distribution1.getValue(i);
					pmf[index][1] = distribution2.getValue(j);
					pmf[index][2] = distribution1.prob(i) * distribution2.prob(j);
					index++;
				}
			return pmf;
		}
		
		return null;
	}
	
	
	
	
}
