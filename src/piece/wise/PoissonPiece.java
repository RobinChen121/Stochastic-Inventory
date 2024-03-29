package piece.wise;

import java.util.Arrays;

import javax.naming.InitialContext;
import javax.swing.JFrame;

import org.apache.commons.math3.ml.distance.DistanceMeasure;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import cern.colt.function.IntIntDoubleFunction;
import umontreal.ssj.probdist.BinomialDist;
import umontreal.ssj.probdist.PoissonDist;

/*
* @author chen
* @date 2022 Nov 7, 22:11:44
* @describe: piecewise approximation for the complementary loss function of the Poisson distribution,
* the heuristic idea is from Roberto Rossi.
*
*/
public class PoissonPiece{
	double lambda;
	int segNum;
	double trunQuantile;
	
	public PoissonPiece(double lambda, int segNum, double trunQuantile) {
		this.lambda = lambda;
		this.segNum = segNum;
		this.trunQuantile = trunQuantile;		
	}
	
	public double[][] partition() {
		PoissonDist dist = new PoissonDist(lambda);
		double[][] result = new double[4][segNum];
		double[] p = new double[segNum];
		for (int i = 0; i < segNum; i++) {
			result[0][i] = 1.0 / (double)segNum; // 
			p[i] = result[0][i];
		}
				
		int supportUB = dist.inverseFInt(trunQuantile);
		int start  = 0;
		int end = 0;
		int index = 0;
		for (int i = 0; i <= supportUB; i++) {
			if (dist.cdf(i) > p[index] + index * p[0] || i == supportUB) {
				end = i == supportUB ? i : i - 1;
				double wOmega = 0;
				for (int j = start; j < end; j++) {
					wOmega += j * dist.prob(j); 
				}
				result[3][index] = wOmega / p[0]; // wOmega is the conditional expectation in a partition Omega, should divide the probilit of this partition supprot
				result[2][index] = end; // wOmega is also the x coordinate of the intersection point of two tangent lines in a partitoin region
				double x = compLossFunction(end) - compLossFunction(start) - end * dist.cdf(end) + start * dist.cdf(start);
				x = x / (dist.cdf(start) - dist.cdf(end));
				result[1][index] = x;
				start = end;
				index++;
			}
		}
		return result;
	}
	
	public void drawPicL(double[] interSecX, double[] endX) {
		PoissonDist dist = new PoissonDist(lambda);
		int supportUB = dist.inverseFInt(trunQuantile);
		int ymax = supportUB;
		int ymin = 0;
		double[] values = new double[ymax - ymin + 1];
		int index = 0;
		int[] support = new int[ymax - ymin + 1];
		for (int y = ymin; y < ymax; y++) {
			values[index] = compLossFunction(y);
			support[index] = y;
			index++;
		}
			
		XYSeries series = new XYSeries("loss function");
		for (int k = 0; k < index; k++) {
			series.add((double) support[k], values[k]);
		}
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
		
//		int start = 0;
//		int end = (int) Math.floor(interSecX[0]);
//		for (int i = 0; i < segNum + 1; i++) {
//			XYSeries series1 = new XYSeries("line");
//			double tangentX = i == 0 ? 0 : endX[i-1]; 		
//			for (int j = start; j < end; j++) {	
//				double y = dist.cdf(tangentX) * (j - tangentX) + compLossFunction((int) tangentX);
//				series1.add((double) j, y);
//			}
//			dataset.addSeries(series1);
//			if (i < segNum) {
//				start = end;
//				end = i == segNum - 1 ? supportUB : (int) Math.floor(interSecX[i+1]);
//			}
//		}
		
		String title = "complementary loss function for Poisson distribution";
		JFreeChart chart = ChartFactory.createXYLineChart(
				title, // chart title
				"y", // x axis label
				"loss function", // y axis label
				dataset, // data
				PlotOrientation.VERTICAL,
				true, // include legend
				false, // tooltips
				false // urls
				);
		XYPlot plot = (XYPlot)chart.getPlot();
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesLinesVisible(0, true); // 设置连线可见
	    plot.setRenderer(renderer);
		
		ChartFrame frame = new ChartFrame("sum of binomial", chart);
		frame.pack(); // fit window to figure size
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	public double compLossFunction(int y) {
		PoissonDist dist = new PoissonDist(lambda);
		double value = 0;
		for (int i = 0; i <= y; i++) {
			value += (y - i) * dist.prob(i);
		}
		return value;
	}
	
	public static void main(String[] args) {
		double lambda = 10;
		int segNum = 4;
		double trunQuantile = 0.999;
		
		PoissonPiece piece = new PoissonPiece(lambda, segNum, trunQuantile);
		double[][] result = piece.partition();
		piece.drawPicL(result[1], result[2]);
		System.out.println(Arrays.deepToString(result));
	}

}
