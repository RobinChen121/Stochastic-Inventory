/**
 * @date: Oct 29, 2020
 */
package milp;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Oct 29, 2020
 * @Desc: codes from Roberto Rossi
 *
 */

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.EmpiricalDist;
import umontreal.ssj.probdist.ExponentialDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.rng.MRG32k3aL;



public class ComplementaryFirstOrderLossFunction {
	Distribution[] distributions;
	long[] seed;
	MRG32k3aL randGenerator;
	
	public ComplementaryFirstOrderLossFunction(Distribution[] distributions, long[] seed){
		this.distributions = distributions;
		this.randGenerator = new MRG32k3aL();
		this.randGenerator.setSeed(seed);
	}
	
	private double[][] sample(int nbSamples){
		this.randGenerator.resetStartStream();
		double[][] sampleMatrix = new double[nbSamples][this.distributions.length];
		for(int i = 0; i < sampleMatrix.length; i++){
			for(int j = 0; j < sampleMatrix[i].length; j++){
				sampleMatrix[i][j] = distributions[j].inverseF(UniformGen.nextDouble(this.randGenerator, 0, 1));
			}
		}
		return sampleMatrix;
	}
	
	
	/**
	 * @param nbSamples
	 * @return the empirical distribution of total demand in the planning horizon
	 * @date: Oct 30, 2020, 5:15:22 PM 
	 */
	public EmpiricalDist getEmpiricalDistribution(int nbSamples){
		double[][] sampleMatrix = this.sample(nbSamples);
		double[] observations = new double[nbSamples]; // row number
		for(int i = 0; i < sampleMatrix.length; i++){ // nbSamples is the row number
			for(int j = 0; j < sampleMatrix[i].length; j++){
				observations[i] += sampleMatrix[i][j];
			}
		}
		Arrays.sort(observations);
		EmpiricalDist empDistribution = new EmpiricalDist(observations);
		return empDistribution;
	}
	
	/**
	 * @param nbSamples
	 * @param precision
	 * @return get the x coordinate and y coordinate for the cdf of the empirical distribution
	 * @date: Oct 30, 2020, 5:21:07 PM 
	 */
	public XYSeries getDistributionXYSeries(int nbSamples, double precision){
		XYSeries series = new XYSeries("Empirical distribution");
		EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples);
		for(int i = 0; i < empDistribution.getN(); i++){
			while(i>0 && i<empDistribution.getN() && empDistribution.getObs(i)==empDistribution.getObs(i-1))i++;
			series.add(empDistribution.getObs(i),empDistribution.cdf(empDistribution.getObs(i)));
		}
		return series;
	}
	
	/**
	 * plot the cdf of the empirical distribution
	 * @param nbSamples
	 * @param precision
	 * @date: Oct 30, 2020, 5:23:23 PM 
	 */
	public void plotEmpiricalDistribution(int nbSamples, double precision){
		XYDataset xyDataset = new XYSeriesCollection(this.getDistributionXYSeries(nbSamples, precision));
		JFreeChart chart = ChartFactory.createXYLineChart("Empirical distribution", "Support", "Frequency",
				 xyDataset, PlotOrientation.VERTICAL, false, true, false);
		ChartFrame frame = new ChartFrame("Empirical distribution",chart);
		frame.setVisible(true);
		frame.setSize(500,400);
	}
	
	public double getLossFunctionValue(double x, int nbSamples){
		EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples);
		double value = 0;
		for(int i = 0; i < empDistribution.getN(); i++){
			value += Math.max(x-empDistribution.getObs(i),0)/empDistribution.getN();
		}
		return value;
	}
	
	public XYSeries getLossFunctionXYSeries(double min, double max, int nbSamples, double precision){
		XYSeries series = new XYSeries("Empirical complementary loss function");
		for(double x = min; x <= max; x+= precision){
			series.add(x, getLossFunctionValue(x, nbSamples));
		}
		return series;
	}
	
	public void plotEmpiricalLossFunction(double min, double max, int nbSamples, double precision, boolean saveToDisk){
		XYSeriesCollection xyDataset = new XYSeriesCollection(this.getLossFunctionXYSeries(min, max, nbSamples, precision));
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
				file.write(lc.toLatex(8, 5));
				file.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public static void main(String[] args){
		testDistributionPlot1();
	}
	
	private static void testDistributionPlot1(){
		long[] seed = {1,2,3,4,5,6};
		Distribution[] distributions = new Distribution[3];
		double lambda[] = {20,5,50};
		distributions[0] = new PoissonDist(lambda[0]);
		distributions[1] = new PoissonDist(lambda[1]);
		distributions[2] = new PoissonDist(lambda[2]);
		ComplementaryFirstOrderLossFunction cfolf = new ComplementaryFirstOrderLossFunction(distributions, seed);
		cfolf.plotEmpiricalDistribution(1000, 1);
	}
	
	private static void testDistributionPlot2(){
		long[] seed = {1,2,3,4,5,6};
		Distribution[] distributions = new Distribution[2];
		distributions[0] = new ExponentialDist(0.1);
		distributions[1] = new NormalDist(10,2);
		ComplementaryFirstOrderLossFunction cfolf = new ComplementaryFirstOrderLossFunction(distributions, seed);
		cfolf.plotEmpiricalDistribution(1000, 1);
	}
	
	private static void testLossFunctionPlot(){
		long[] seed = {1,2,3,4,5,6};
		Distribution[] distributions = new Distribution[3];
		double lambda[] = {20,5,50};
		distributions[0] = new PoissonDist(lambda[0]);
		distributions[1] = new PoissonDist(lambda[1]);
		distributions[2] = new PoissonDist(lambda[2]);
		ComplementaryFirstOrderLossFunction cfolf = new ComplementaryFirstOrderLossFunction(distributions, seed);
		cfolf.plotEmpiricalLossFunction(50, 90, 1000, 1, true);
	}
}
