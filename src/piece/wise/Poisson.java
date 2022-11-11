package piece.wise;

import umontreal.ssj.probdist.PoissonDist;

/*
* @author chen
* @date 2022 Nov 7, 22:11:44
* @describe: piecewise approximation for the complementary loss function of the Poisson distribution
*
*/
public class Poisson {
	public double[][] partition(double lambda, int segNum) {
		PoissonDist dist = new PoissonDist(lambda);
		double[][] result = new double[2][segNum];
		for (int i = 0; i < segNum; i++)
			result[0][i] = 1.0 / (double)segNum;
		return null;
	}
	
	public static void main(String[] args) {
		double lambda = 10;
		int segNum = 4;
	}

}
