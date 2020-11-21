/**
 * @date: Oct 31, 2020
 */
package milp;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Oct 31, 2020
 * @Desc: codes are from Roberto Rossi
 *
 */

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.EmpiricalDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;



public class PiecewiseComplementaryFirstOrderLossFunction extends
ComplementaryFirstOrderLossFunction {

	public PiecewiseComplementaryFirstOrderLossFunction(Distribution[] distributions, long[] seed){
		super(distributions, seed);
	}
	
	
	
	/**
	 * @param probabilityMasses
	 * @param nbSamples
	 * @return get conditional expectations for the empirical distributions from the samples of the distributions,
	 * based on the definition of conditional expectations for discrete distribution
	 * @date: Nov 1, 2020, 1:49:13 PM 
	 */
	public double[] getConditionalExpectations(double[] probabilityMasses, int nbSamples){
		double[] conditionalExpectations = new double[probabilityMasses.length];
		EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples);
		double probabilityMass = 0;
		int conditionalExpectationIndex = 0;
		// getN(): get the number of observations which is also the number of samples
		// getObs(i): get the value of X_i
		int n =  empDistribution.getN();
		for(int i = 0; i < empDistribution.getN(); i++){
			if(probabilityMass < 1 && probabilityMass < probabilityMasses[conditionalExpectationIndex]){
				conditionalExpectations[conditionalExpectationIndex] += empDistribution.getObs(i)/empDistribution.getN();
				probabilityMass += 1.0/empDistribution.getN();
			}else{
				conditionalExpectations[conditionalExpectationIndex] /= probabilityMasses[conditionalExpectationIndex];
				probabilityMass = 0;
				conditionalExpectationIndex++;
			}
		}
		conditionalExpectations[conditionalExpectationIndex] /= probabilityMasses[conditionalExpectationIndex];
		return conditionalExpectations;
	}

	public XYSeries getLossFunctionXYSeriesForSegment(int segmentIndex, 
			double[] probabilityMasses, 
			double[] conditionalExpectations, 
			double min, 
			double max, 
			double minYValue,
			double precision){
		XYSeries series = new XYSeries("Piecewise complementary loss function" + Integer.toString(segmentIndex));
		for(double x = min; x <= max; x+= precision){
			double value = getPiecewiseLossFunctionValue(segmentIndex, x, probabilityMasses, conditionalExpectations);
			if(value >= minYValue) series.add(x, value);
		}
		return series;
	}

	/**
	 * @param segmentIndex
	 * @param x
	 * @param probabilityMasses
	 * @param conditionalExpectations
	 * @return get the piecewise value of the (complementary) loss function (I^+)
	 * @date: Nov 1, 2020, 4:56:50 PM 
	 */
	public double getPiecewiseLossFunctionValue(int segmentIndex, double x, double[] probabilityMasses, double[] conditionalExpectations){
		double value = 0;
		for(int i = 1; i <= segmentIndex; i++){
			value += (x-conditionalExpectations[i-1])*probabilityMasses[i-1];
		}
		return value;
	}
	
	public XYSeries getPiecewiseErrorXYSeries(double min, double max, int nbSamples, double[] probabilityMasses, double[] conditionalExpectations, double precision){
		XYSeries series = new XYSeries("Piecewise complementary loss function error");
		for(double x = min; x <= max; x+= precision){
			series.add(x, getPiecewiseErrorValue(x, nbSamples, probabilityMasses, conditionalExpectations));
		}
		return series;
	}
	
	public double getPiecewiseErrorValue(double x, int nbSamples, double[] probabilityMasses, double[] conditionalExpectations){
		double lossFunctionValue = this.getLossFunctionValue(x, nbSamples); // loss function value of this x

		double maxValue = 0;
		// choose a maximum error value for all the segments
		for(int j = 0; j <= probabilityMasses.length; j++){
			double value = 0;
			for(int i = 1; i <= j; i++){
				value += (x-conditionalExpectations[i-1])*probabilityMasses[i-1];
			}
			maxValue = Math.max(maxValue, value);
		}
		
		return lossFunctionValue-maxValue;
	}

	
	/**
	 * draw the piecewise loss function
	 * @param min
	 * @param max
	 * @param minYValue
	 * @param probabilityMasses
	 * @param nbSamples
	 * @param precision
	 * @param saveToDisk
	 * @date: Oct 31, 2020, 5:18:12 PM 
	 */
	public void plotPiecewiseLossFunction(double min, double max, double minYValue, double[] probabilityMasses, int nbSamples, double precision, boolean saveToDisk){
		int segments = probabilityMasses.length + 1;
		double[] conditionalExpectations = this.getConditionalExpectations(probabilityMasses, nbSamples);

		XYSeriesCollection xyDataset = new XYSeriesCollection();

		xyDataset.addSeries(this.getLossFunctionXYSeries(min, max, nbSamples, precision));
		
		// draw each segment
		for(int i = 0; i < segments; i++)
			xyDataset.addSeries(this.getLossFunctionXYSeriesForSegment(i, probabilityMasses, conditionalExpectations, min, max, minYValue, precision));

		xyDataset.addSeries(this.getPiecewiseErrorXYSeries(min, max, nbSamples, probabilityMasses, conditionalExpectations, precision));
		
		JFreeChart chart = ChartFactory.createXYLineChart("Empirical complementary loss function", "x", "CL(x)",
				xyDataset, PlotOrientation.VERTICAL, false, true, false);
		ChartFrame frame = new ChartFrame("Empirical complementary loss function",chart);
		frame.setVisible(true);
		frame.setSize(500,400);

		if(saveToDisk){

			XYLineChart lc = new XYLineChart("Piecewise linearization", "x", "CL(x)", xyDataset);

			try {
				File latexFolder = new File("./latex");
				if(!latexFolder.exists()){
					latexFolder.mkdir();
				}
				Writer file = new FileWriter("./latex/graph.tex");
				file.write(lc.toLatex(5, 5));
				file.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private static double[] toPrimitive(Double[] array) {
		   if (array == null) {
		     return null;
		   } else if (array.length == 0) {
		     return new double[0];
		   }
		   final double[] result = new double[array.length];
		   for (int i = 0; i < array.length; i++) {
		     result[i] = array[i].doubleValue();
		   }
		   return result;
		 }

	/**
	 * @param probabilityMasses
	 * @param nbSamples
	 * @return approximation errors for each segment
	 * @date: Nov 1, 2020, 7:21:26 PM 
	 */
	public double[] getApproximationErrors(double[] probabilityMasses, int nbSamples){
		double[] conditionalExpectations = this.getConditionalExpectations(probabilityMasses, nbSamples);
		double[] approximationErrors = new double[conditionalExpectations.length];
		
		// max error locates at break points where the x coordinate is conditionalExpectations[i]
		for(int i = 0; i < probabilityMasses.length; i++){
			approximationErrors[i] = getLossFunctionValue(conditionalExpectations[i], nbSamples)-
					getPiecewiseLossFunctionValue(i, conditionalExpectations[i], probabilityMasses, conditionalExpectations);
		}

		return approximationErrors;
	}

	public double getMaxApproximationError(double[] probabilityMasses, int nbSamples){
		double[] approximationErrors = this.getApproximationErrors(probabilityMasses, nbSamples);
		double maxApproximationError = 0;

		for(int i = 0; i < probabilityMasses.length; i++){
			maxApproximationError = Math.max(maxApproximationError, approximationErrors[i]);
		}

		return maxApproximationError;
	}

	private static void testApproximationErrors(){
		long[] seed = {1,2,3,4,5,6};
		Distribution[] distributions = new Distribution[1];
		distributions[0] = new PoissonDist(20);
		//distributions[1] = new ExponentialDist(0.1);
		PiecewiseComplementaryFirstOrderLossFunction pwcfolf = new PiecewiseComplementaryFirstOrderLossFunction(distributions, seed);
		double[] probabilityMasses = new double[10];
		Arrays.fill(probabilityMasses, 0.1); // assume all p_i are equal
		System.out.println("probability p_i for each segment: ");
		System.out.println(Arrays.toString(probabilityMasses));
		int nbSamples = 1000;
		double[] approximationErrors = pwcfolf.getApproximationErrors(probabilityMasses, nbSamples);
		System.out.println("conditional expectations in each segement: ");
		double[] conditionExpects = pwcfolf.getConditionalExpectations(probabilityMasses, nbSamples);
		System.out.println(Arrays.toString(conditionExpects));
		System.out.println("approximation error for each segment: ");
		for(int i = 0; i < probabilityMasses.length; i++){
			System.out.print(approximationErrors[i]+"\t");
		}
	}
	
	private static void testPiecewiseLossFunction(){
		long[] seed = {1,2,3,4,5,6};
		Distribution[] distributions = new Distribution[1];
		//distributions[0] = new NormalDist(0,1);
		//distributions[1] = new ExponentialDist(0.1);
		distributions[0] = new PoissonDist(20);
		PiecewiseComplementaryFirstOrderLossFunction pwcfolf = new PiecewiseComplementaryFirstOrderLossFunction(distributions, seed);
		double[] probabilityMasses = {0.2, 0.2, 0.2, 0.2, 0.2};
		int nbSamples = 1000;
		pwcfolf.plotPiecewiseLossFunction(-2, 2, -1, probabilityMasses, nbSamples, 0.1, false);
	}
	


	public static void main(String[] args){
		//testPiecewiseLossFunction();
		testApproximationErrors();
	}
}
