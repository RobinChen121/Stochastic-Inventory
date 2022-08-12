package test;

import static org.junit.Assert.assertArrayEquals;

import org.apache.commons.math3.util.Precision;
import org.junit.Test;

import nonlinear.optimization.Calcfc;
import nonlinear.optimization.Cobyla;
import nonlinear.optimization.CobylaExitStatus;

import java.util.Arrays;

public class TestJcobyla {
	// FIELDS
    
    private double rhobeg = 0.5;
    private double rhoend = 1.0e-6;
    private int iprint = 1;
    private int maxfun = 3500;

    // TESTS
    
    /**
     * Minimization of a simple quadratic function of two variables.
     */
    public void test01FindMinimum() {
        System.out.format("%nOutput from test problem 1 (Simple quadratic)%n");
        Calcfc calcfc = new Calcfc() {
            public double compute(int n, int m, double[] x, double[] con) {
                return 10.0 * Math.pow(x[0] + 1.0, 2.0) + Math.pow(x[1], 2.0);
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.findMinimum2(calcfc, 2, 0, x, rhobeg, rhoend, iprint, maxfun);
        System.out.println(result);
        assertArrayEquals(null, new double[] { -1.0, 0.0 }, x, 1.0e-5);
    }
    
    /**
     * Easy two dimensional minimization in unit circle.
     */
    @Test
    public void test02FindMinimum() {
        System.out.format("%nOutput from test problem 2 (2D unit circle calculation)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
                con[0] = 1.0 - x[0] * x[0] - x[1] * x[1]; // >=
                con[1] = -x[0] * x[0] - x[1] * x[1] + 1.0;
                return x[0] * x[1]; // objective
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.findMinimum2(calcfc, 2, 2, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { Math.sqrt(0.5), -Math.sqrt(0.5) }, x, 1.0e-5);
    }

    /**
     * Easy three dimensional minimization in ellipsoid.
     */
    @Test
    public void test03FindMinimum() {
        System.out.format("%nOutput from test problem 3 (3D ellipsoid calculation)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
                con[0] = 1.0 - x[0] * x[0] - 2.0 * x[1] * x[1] - 3.0 * x[2] * x[2];
                return x[0] * x[1] * x[2];
            }
        };
        double[] x = {1.0, 1.0, 1.0 };
        CobylaExitStatus result = Cobyla.findMinimum2(calcfc, 3, 1, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { 1.0 / Math.sqrt(3.0), 1.0 / Math.sqrt(6.0), -1.0 / 3.0 }, x, 1.0e-5);
    }
    
    /**
     * This problem is taken from Fletcher's book Practical Methods of
     * Optimization and has the equation number (9.1.15).
     */
    @Test
    public void test06FindMinimum() {
        System.out.format("%nOutput from test problem 6 (Equation (9.1.15) in Fletcher's book)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
                con[0] = x[1] - x[0] * x[0]; // x_2 - x_1^2 >= 0
                con[1] = 1.0 - x[0] * x[0] - x[1] * x[1];
                return -x[0] - x[1];
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.findMinimum2(calcfc, 2, 2, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { Math.sqrt(0.5), Math.sqrt(0.5) }, x, 1.0e-5);
    }
    
    /**
     * This problem is taken from Fletcher's book Practical Methods of
     * Optimization and has the equation number (14.4.2).
     */
    @Test
    public void test07FindMinimum() {
        System.out.format("%nOutput from test problem 7 (Equation (14.4.2) in Fletcher)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
                con[0] = 5.0 * x[0] - x[1] + x[2]; // z + 5x_1 - x_2 >= 0
                con[1] = x[2] - x[0] * x[0] - x[1] * x[1] - 4.0 * x[1];
                con[2] = x[2] - 5.0 * x[0] - x[1];
                return x[2];
            }
        };
        double[] x = {1.0, 1.0, 1.0 };
        double[] result = Cobyla.findMinimum(calcfc, 3, 3, x, rhobeg, rhoend, iprint, maxfun);
        double[] result2 = Arrays.stream(result).map(i -> Precision.round(i, 1)).toArray();
        System.out.println("X = " + Arrays.toString(result2));
        assertArrayEquals(null, new double[] { 0.0, -3.0, -3.0 }, x, 1.0e-5);
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

    public void testScenario() {
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
                	con[i + 2] = pvalues[i];
                }
                
                double mean = mean(xvalues);
                double variance = variance(xvalues, pvalues);                             
                return Math.pow(mean - 5, 2) + Math.pow(Math.sqrt(variance) - 2, 2);
            }
        };
        double[] x = {2.0, 5.0, 5.0, 8.0, 5.0, 0.2, 0.2, 0.2, 0.2, 0.2};
        int halfLength = x.length/2;
        double[] result = Cobyla.findMinimum(calcfc, x.length, halfLength + 2, x, rhobeg, rhoend, iprint, maxfun);
        double[] result2 = Arrays.stream(result).map(i -> Precision.round(i, 2)).toArray();
        System.out.println("X = " + Arrays.toString(result2));
        
    	double[] pvalues = new double[halfLength];
    	double[] xvalues = new double[halfLength];
    	System.arraycopy(result2, 0, xvalues, 0, halfLength);
    	System.arraycopy(result2, halfLength, pvalues, 0, halfLength);
    	
        System.out.println("mean = " + mean(xvalues));
        System.out.println("variance = " + variance(xvalues, pvalues));
    }
    
	public static void main(String[] args) {
		TestJcobyla Test = new TestJcobyla();
		Test.testScenario();

	}

}
