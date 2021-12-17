package sdp.chance;

import java.util.Arrays;
import java.util.Iterator;

import umontreal.ssj.probdist.BinomialDist;

public class ComputeStatisticUpperBound {
	
	static double getUpperBound(double epsilon, double tau, int M, int N, double[] probs) {
		int p = (int) (N*epsilon);
		double theta = 0;
		for (int i = 0; i < p; i++) {
			theta += BinomialDist.prob(N, epsilon, i);
		}
		int L = 0;
		Arrays.sort(probs);
		for (int i = 0; i < M; i++) {
			double pi = BinomialDist.prob(M, theta, i);
			if (pi > tau) {
				L = i - 1;
				break;
			}
		}
		return probs[L];
	}
	

	public static void main(String[] args) {
		double[] probs = new double[] {0.7659, 0.7500, 0.8323, 0.8329, 0.8124, 0.7526, 0.7570, 0.8954};
		int N = 625;
		double epsilon = 0.4; 
		double tau = 0.05;
		int M = probs.length;
		double a = getUpperBound(epsilon, tau, M, N, probs);
		System.out.println(a);
		

	}

}
