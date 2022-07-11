package workforce;

import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;


public class GetPmf {
	double truncationQuantile;
	
	public GetPmf(double truncationQuantile) {
		this.truncationQuantile = truncationQuantile;
	}
	
	public double[][] getPmf(Distribution distribution){
		double supportLB;
		double supportUB;
		
		supportLB = 0;
		supportUB = (int) distribution.inverseF(truncationQuantile);
		int demandLength = (int) (supportUB - supportLB + 1);
		double[][] pmf = new double[demandLength][];
		
		for (int j = 0; j < demandLength; j++) {
			pmf[j][0] = supportLB + j;		
			double probilitySum = distribution.cdf(supportUB) - distribution.cdf(supportLB);
			pmf[j][1] = ((DiscreteDistributionInt) distribution).prob(j) / probilitySum;
		}
		return pmf;
	}

}
