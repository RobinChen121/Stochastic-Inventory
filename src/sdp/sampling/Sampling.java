package sdp.sampling;

import java.util.Arrays;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.probdistmulti.BiNormalDist;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.randvar.UniformIntGen;
import umontreal.ssj.rng.MRG31k3p;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

/** 
* @author chen zhen 
* @version 2018, April 11th, 4:18:55 pm 
* @value random sampling and latin hypercube sampling;
* 
* give same seeds if you want to get the same random number in different running time,
* it seems the seeds in the ssj are complex which needs 6 seeds to generate random number;
*/

public class Sampling {
	
	static RandomStream stream = new MRG32k3a();
	
	/**
	 * Reinitializes the stream to its initial state.
	 */
	public static void resetStartStream(){
		stream.resetStartStream();
	}
	
	
	/**
	 * Reinitializes the stream to the beginning of its next substream.
	 */
	public static void resetNextSubstream(){
		stream.resetNextSubstream();
	}
	
	/** random sampling
	 * @param distributions
	 * @param sampleNum
	 * @return a 2D random samples
	 */
	public double[][] generateRanSamples(Distribution[] distributions, int sampleNum){
		int periodNum = distributions.length;		
		double[][] samples = new double[sampleNum][periodNum]; 
		
		for (int i = 0; i < periodNum; i++)
			for (int j = 0; j < sampleNum; j++) {
				samples[j][i] =  UniformGen.nextDouble(stream, 0, 1.0);
				samples[j][i] = distributions[i].inverseF(samples[j][i]);
			}
		return samples;
	}
	
	/**
	 * 
	 * @param distributions
	 * 
	 * @return one demands sample
	 */
	
	public double[] getNextSample(Distribution[] distributions) {
		int periodNum = distributions.length;
		double[] sample = new double[periodNum];
		UniformGen uniform = new UniformGen(stream);
		for (int i = 0; i < periodNum; i++) {
			sample[i] = distributions[i].inverseF(uniform.nextDouble());
		}
		return sample;
	}
	
	
	
	/** latin hypercube sampling
	 * @param distributions
	 * @param sampleNum
	 * @return a 2D random samples
	 */
	public double[][] generateLHSamples(Distribution[] distributions, int sampleNum){
		int periodNum = distributions.length;		
		double[][] samples = new double[sampleNum][periodNum]; 
		
		// generate random possibility in [i/n, (i+1)/n], then get percent point function according to the possibility		
		for (int i = 0; i < periodNum; i++)
			for (int j = 0; j < sampleNum; j++) {
				double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
				double lowBound = (double) j/ (double) sampleNum;
				samples[j][i] = lowBound + randomNum;
				samples[j][i] = distributions[i].inverseF(samples[j][i]);
			}
		
	    shuffle(samples); // ��������
		return samples;
	}
	
	/** latin hypercube sampling
	 * @param distributions
	 * @param sampleNum
	 * @return a 2D random samples, in which the sample number in each period can be different;
	 * each row is a period
	 */
	public double[][] generateLHSamples(Distribution[] distributions, int[] sampleNums){
		resetNextSubstream();
		
		int periodNum = distributions.length;		
		double[][] samples = new double[periodNum][]; 
		
		// generate random possibility in [i/n, (i+1)/n], then get percent point function according to the possibility		
		for (int i = 0; i < periodNum; i++) {
			int sampleNum = sampleNums[i];
			samples[i] = new double[sampleNum];
			for (int j = 0; j < sampleNum; j++) {
				double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
				double lowBound = (double) j/ (double) sampleNum;
				samples[i][j] = lowBound + randomNum;
				samples[i][j] = distributions[i].inverseF(samples[i][j]);
			}
		}
		
	    shuffle2(samples); // ��������
		return samples;
	}
	
	/** latin hypercube sampling for binormal distribution.
	 * 
	 * Since two independent variable, generate two variable independently, and merge the two samples into one
	 * @param distributions
	 * @param sampleNum
	 * @return a 2D random samples
	 */
	public double[][] generateLHSamples(BiNormalDist[] distributions, int sampleNum){
		int periodNum = distributions.length;		
		double[][] samples = new double[sampleNum][periodNum * 2]; 
		
		double[][] samples1 = new double[sampleNum][periodNum];
		double[][] samples2 = new double[sampleNum][periodNum];
		
		// generate random possibility in [i/n, (i+1)/n], then get percent point function according to the possibility		
		for (int i = 0; i < periodNum; i++) {
			NormalDist distribution1 = new NormalDist(distributions[i].getMu1(), distributions[i].getSigma1());
			for (int j = 0; j < sampleNum; j++) {
				double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
				double lowBound = (double) j/ (double) sampleNum;
				samples1[j][i] = lowBound + randomNum;
				samples1[j][i] = distribution1.inverseF(samples1[j][i]);
			}		
			shuffle(samples1); // ��������		
		}
		for (int i = 0; i < periodNum; i++) {
			NormalDist distribution2 = new NormalDist(distributions[i].getMu2(), distributions[i].getSigma2());
			for (int j = 0; j < sampleNum; j++) {
				double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
				double lowBound = (double) j/ (double) sampleNum;
				samples2[j][i] = lowBound + randomNum;
				samples2[j][i] = distribution2.inverseF(samples2[j][i]);
			}		
			shuffle(samples2); // ��������		
		}
		
		for (int i = 0; i < sampleNum; i++) {
			for (int j = 0; j < periodNum; j++) {
				samples[i][j] = samples1[i][j];
				samples[i][j + periodNum] = samples2[i][j];
			}			
		}
		
		return samples;
	}
	
	
//	/** latin hypercube sampling for bi poisson distribution.
//	 * 
//	 * Since two independent variable, generate two variable independently, and merge the two samples into one
//	 * @param distributions
//	 * @param sampleNum
//	 * @return a 2D random samples
//	 */
//	public double[][] generateLHSamplesMulti(Distribution[][] distributions, int sampleNum){
//		int periodNum = distributions.length;		
//		double[][] samples = new double[sampleNum][periodNum * 2]; 
//		
//		double[][] samples1 = new double[sampleNum][periodNum];
//		double[][] samples2 = new double[sampleNum][periodNum];
//		
//		// generate random possibility in [i/n, (i+1)/n], then get percent point function according to the possibility		
//		for (int i = 0; i < periodNum; i++) {
//			Distribution distribution1 = distributions[i][0];
//			for (int j = 0; j < sampleNum; j++) {
//				double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
//				double lowBound = (double) j/ (double) sampleNum;
//				samples1[j][i] = lowBound + randomNum;
//				samples1[j][i] = distribution1.inverseF(samples1[j][i]);
//			}		
//			shuffle(samples1); // ��������		
//		}
//		for (int i = 0; i < periodNum; i++) {
//			Distribution distribution2 = distributions[i][1];
//			for (int j = 0; j < sampleNum; j++) {
//				double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
//				double lowBound = (double) j/ (double) sampleNum;
//				samples2[j][i] = lowBound + randomNum;
//				samples2[j][i] = distribution2.inverseF(samples2[j][i]);
//			}		
//			shuffle(samples2); // ��������		
//		}
//		
//		for (int i = 0; i < sampleNum; i++) {
//			for (int j = 0; j < periodNum; j++) {
//				samples[i][j] = samples1[i][j];
//				samples[i][j + periodNum] = samples2[i][j];
//			}			
//		}
//		
//		return samples;
//	}
//	
	/** latin hypercube sampling with truncationQuantile 
	 * @param distributions
	 * @param sampleNum
	 * @return a 2D random samples
	 */
	public double[][] generateLHSamples(Distribution[] distributions, int sampleNum, double frac){
		int periodNum = distributions.length;		
		double[][] samples = new double[sampleNum][periodNum]; 
		
		// ��ÿ��[i/n, (i+1)/n] ������һ��������ʣ�Ȼ����ݸ��ʵõ�ָ���ֲ�����
		for (int i = 0; i < periodNum; i++)
			for (int j = 0; j < sampleNum; j++) {
				double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
				double lowBound = (double) j/ (double) sampleNum;
				samples[j][i] = frac*(lowBound + randomNum);
				samples[j][i] = distributions[i].inverseF(samples[j][i]);
			}
		
	    shuffle(samples); // ��������
		return samples;
	}
	
	/** shuffle a 2D array
	 * 
	 */
	 double[][] shuffle(double[][] samples){
		for(int i = 0; i < samples[0].length; i++)
			for (int j = 0; j < samples.length; j++){
			int mark = UniformIntGen.nextInt(stream, 0, samples.length - 1);
			double temp = samples[j][i];
			samples[j][i] = samples[mark][i];
			samples[mark][i] = temp;
		}
		return samples;
	}
	 
	 /** shuffle a 2D array
		* 
	 */
	 double[][] shuffle2(double[][] samples){
		for(int i = 0; i < samples.length; i++)  // t
			for (int j = 0; j < samples[i].length; j++){
			int mark = UniformIntGen.nextInt(stream, 0, samples[i].length - 1);
			double temp = samples[i][j];
			samples[i][j] = samples[i][mark];
			samples[i][mark] = temp;
		}
		return samples;
	}
	 
		/** latin hypercube sampling
		 * @param distributions
		 * @param sampleNum
		 * @param t
		 * @return a 2D random samples for multi products in a period t, very slow for multi products
		 * @date: Apr 29, 2020, 12:04:34 PM 
		 */
		public double[][] generateLHSamplesMultiProd(Distribution[][] distributions, int sampleNum){
			int itemNum = distributions.length;		
			int T = distributions[0].length;
			double[][] samples = new double[sampleNum][itemNum * T]; 
			
			// generate random possibility in [i/n, (i+1)/n], then get percent point function according to the possibility		
			for (int i = 0; i < sampleNum; i++) {
				for (int j = 0; j < itemNum; j++)
					for (int t = 0; t < T; t++) {
						double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
						double lowBound = (double) i / (double) sampleNum;
						double ppf = lowBound + randomNum;
						samples[i][j + t * itemNum] = Math.round(distributions[j][t].inverseF(ppf) * 1.0) / 1.0;
					}
			}
		    shuffle(samples); // ��������
			return samples;
		}
	
	
//	public static void main(String[] args) {
//		int sampleNum = 1000;
//		double[] meanDemand = {20, 5};	
//		
//		PoissonDist[] distributions = new PoissonDist[meanDemand.length];
//		for (int i = 0; i < meanDemand.length; i++)
//			distributions[i] = new PoissonDist(meanDemand[i]);
//		
//		double[][] samples = Sampling.generateRanSamples(distributions, sampleNum);
//		System.out.println(Arrays.deepToString(samples));
//		
//		for (int i = 0; i < meanDemand.length; i++) {
//			double sum = 0; 
//			for (int j = 0; j < sampleNum; j++) {
//				sum += samples[j][i];
//			}
//			System.out.println(sum/sampleNum);
//		}	
//	}
}
