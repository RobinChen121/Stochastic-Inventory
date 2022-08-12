package nonlinear.optimization;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.DoubleStream;

import javax.annotation.processing.Generated;

import org.apache.commons.math3.util.Precision;

import cern.colt.function.IntIntDoubleFunction;
import milp.GenerateScenarioTree;

/**
 * @author chen
 * @date ：2022 July 25, 22:38:27 
 * @ the nonlinear optimization package in java is not very suitable to generate scenarios, probabilities so close
 *
 */
public class GenearateScenarios {
	// FIELDS  
    private double rhobeg = 0.5;
    private double rhoend = 1.0e-6;
    private int iprint = 0;
    private int maxfun = 3500;
    
	double[] mean;
	double[] std;
	int sampleNum;
	
	public GenearateScenarios(double[] mean, double[] std, int sampleNum) {
		this.mean = mean;
		this.std = std;
		this.sampleNum = sampleNum;	
	}
	
	double mean(double[] xvalues) {
    	return Arrays.stream(xvalues).average().getAsDouble();
    }
    
    double variance(double[] xvalues, double[] pvalues) {
        double variance = 0;
        double mean = mean(xvalues);
        for(int i = 0; i < xvalues.length; i++) {
        	variance += pvalues[i] * Math.pow(xvalues[i] - mean, 2);
        }
        return variance;
    }
    
	public void generateInOnePeriod(double giveMean, double giveStd) {
        Calcfc calcfc = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
            	int halfLength = x.length/2;
            	double[] pvalues = new double[halfLength];
            	double[] xvalues = new double[halfLength];
            	System.arraycopy(x, 0, xvalues, 0, halfLength);
            	System.arraycopy(x, halfLength, pvalues, 0, halfLength);
            	
                con[0] = Arrays.stream(pvalues).sum() - 1; 
                con[1] = -Arrays.stream(pvalues).sum() + 1; 
                
                for(int i = 0; i < xvalues.length; i++) {
                	con[i + 2] = pvalues[i] - 0.1;
                }
                
                double mean = mean(xvalues);
                double variance = variance(xvalues, pvalues);                             
                return Math.pow(mean - giveMean, 2) + Math.pow(variance - Math.pow(giveStd, 2), 2);
            }
        };
        double iniProb = 1.0 / (double) sampleNum;
        double[] x = new double[2*sampleNum];
        for (int i = 0; i < sampleNum; i++){
        	x[sampleNum + i] = iniProb;
        	x[i] = giveMean + Math.pow(-1, i) * Math.random()*giveStd;
        }
        int halfLength = x.length/2;
        double[] result = Cobyla.findMinimum(calcfc, x.length, halfLength + 2, x, rhobeg, rhoend, iprint, maxfun);
        double[] result2 = Arrays.stream(result).map(i -> Precision.round(i, 2)).toArray();
        System.out.println("X = " + Arrays.toString(result2));
        
    	double[] pvalues = new double[halfLength];
    	double[] xvalues = new double[halfLength];
    	System.arraycopy(result2, 0, xvalues, 0, halfLength);
    	System.arraycopy(result2, halfLength, pvalues, 0, halfLength);
    	
        System.out.println("mean = " + mean(xvalues));
        System.out.println("std = " + Math.sqrt(variance(xvalues, pvalues)));
	}
	
	
	public double[][][] generateTree(){
		int n = mean.length;
		double[][][] result = new double[2][n][sampleNum];
		for (int i = 0; i < n; i++) {
			
		}
		return null;
	}

	public static void main(String[] args) {
		double[] meanDemand = {63, 27, 10, 24, 18};
		double[] sigma = DoubleStream.of(meanDemand).map(i -> i*0.25).toArray();
		 GenearateScenarios test = new GenearateScenarios(meanDemand, sigma, 3);
		 test.generateInOnePeriod(meanDemand[0], sigma[0]);
		 

	}

}
