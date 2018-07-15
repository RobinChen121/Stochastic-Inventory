package sdp.sampling;

import java.util.Arrays;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.randvar.UniformIntGen;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

/** 
* @author chen zhen 
* @version 创建时间：2018年4月11日 下午4:18:55 
* @value 类说明: random sampling and latin hypercube sampling
*/

public class Sampling {
	
	static RandomStream stream = new MRG32k3aL();
	
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
	public static double[][] generateRanSamples(Distribution[] distributions, int sampleNum){
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
	 * @return one demands sample
	 */
	
	public static double[] getNextSample(Distribution[] distributions) {
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
	public static double[][] generateLHSamples(Distribution[] distributions, int sampleNum){
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
		
	    shuffle(samples); // 打乱数组
		return samples;
	}
	
	/** latin hypercube sampling with truncationQuantile 
	 * @param distributions
	 * @param sampleNum
	 * @return a 2D random samples
	 */
	public static double[][] generateLHSamples(Distribution[] distributions, int sampleNum, double frac){
		int periodNum = distributions.length;		
		double[][] samples = new double[sampleNum][periodNum]; 
		
		// 在每个[i/n, (i+1)/n] 内生成一个随机概率，然后根据概率得到指定分布的数
		for (int i = 0; i < periodNum; i++)
			for (int j = 0; j < sampleNum; j++) {
				double randomNum = UniformGen.nextDouble(stream, 0, 1.0/sampleNum);
				double lowBound = (double) j/ (double) sampleNum;
				samples[j][i] = frac*(lowBound + randomNum);
				samples[j][i] = distributions[i].inverseF(samples[j][i]);
			}
		
	    shuffle(samples); // 打乱数组
		return samples;
	}
	
	/** shuffle a 2D array
	 * 
	 */
	static double[][] shuffle(double[][] samples){
		for(int i = 0; i < samples[0].length; i++)
			for (int j = 0; j < samples.length; j++){
			int mark = UniformIntGen.nextInt(stream, 0, samples.length - 1);
			double temp = samples[j][i];
			samples[j][i] = samples[mark][i];
			samples[mark][i] = temp;
		}
		return samples;
	}
	
	
	public static void main(String[] args) {
		int sampleNum = 1000;
		double[] meanDemand = {20, 5};	
		
		PoissonDist[] distributions = new PoissonDist[meanDemand.length];
		for (int i = 0; i < meanDemand.length; i++)
			distributions[i] = new PoissonDist(meanDemand[i]);
		
		double[][] samples = Sampling.generateRanSamples(distributions, sampleNum);
		System.out.println(Arrays.deepToString(samples));
		
		for (int i = 0; i < meanDemand.length; i++) {
			double sum = 0; 
			for (int j = 0; j < sampleNum; j++) {
				sum += samples[j][i];
			}
			System.out.println(sum/sampleNum);
		}	
	}
}
