package sdp.sampling;

import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.randvar.UniformIntGen;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

public class SampleFactory {
	
	private static RandomStream stream = new MRG32k3aL();
	
	/**
	 * Reinitializes the stream to the beginning of its next substream.
	 */
	public void resetNextSubstream(){
		stream.resetNextSubstream();
	}
	
	/**
	 * Reinitializes the stream to its initial state.
	 */
	public void resetStartStream(){
		stream.resetStartStream();
	}
	
	/**
	 * Implements Simple Random Sampling
	 * @param distributions array of distributions to be sampled
	 * @return a Simple Random Sample for the distributions in {@code distributions}
	 */
	public static double[] getNextSample(Distribution[] distributions){
		UniformGen uniform = new UniformGen(stream);
		return IntStream.iterate(0, i -> i + 1)
		                .limit(distributions.length).mapToDouble(
		                      i -> distributions[i].inverseF(uniform.nextDouble()))
		                .toArray();
	}
	
	/**
	 * Implements Latin Hypercube Sampling as originally introduced in 
	 * 
	 * McKay, M.D.; Beckman, R.J.; Conover, W.J. (May 1979). 
	 * "A Comparison of Three Methods for Selecting Values of Input Variables 
	 * in the Analysis of Output from a Computer Code". 
	 * Technometrics 21 (2): 
	 * 
	 * @param distributions array of distributions to be sampled 
	 * @param samples number of samples
	 * @return a Latin Hypercube Sample for the distributions in {@code distributions}
	 */
	public double[][] getNextLHSample(Distribution[] distributions, int samples){
		double x[][] = new double[distributions.length][samples];
		x = IntStream.iterate(0, d -> d + 1)
		             .limit(distributions.length)
		             .mapToObj(
		                   d -> DoubleStream.iterate(0, i -> i + 1.0/samples)
		                                    .limit(samples)
		                                    .map(i -> distributions[d].inverseF(i + UniformGen.nextDouble(stream, 0, 1.0/samples)))
		                                    .toArray())
		             .toArray(double[][]::new);	
		for(int i = 0; i < x.length; i++){
			shuffle(x[i]);
		}
		return x;
	}
	
	/**
	 * Returns a random shuffle of {@code sample}.
	 * 
	 * @param sample the original sample.
	 * @return a random shuffle of {@code sample}.
	 */
	private static double[] shuffle(double[] sample){
		for(int i = 0; i < sample.length; i++){
			int j = UniformIntGen.nextInt(stream, 0, sample.length - 1);
			double temp = sample[i];
			sample[i] = sample[j];
			sample[j] = temp;
		}
		return sample;
	}
}
