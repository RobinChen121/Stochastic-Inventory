package piece.wise;

import java.util.Arrays;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import umontreal.ssj.probdist.BinomialDist;

/*
* @author chen
* @date 2022 Dec 1, 13:21:23
* @describe: piece wise for binomial distribution
*
*/
public class BinomialPiece {
	int ymin;
	int ymax;
	int w;
	double p;
	int segNum;
	
	public BinomialPiece(int ymin, int ymax, int w, double p, int segNum){
		this.ymin = ymin;
		this.ymax = ymax;
		this.w = w;
		this.p = p;		
		this.segNum = segNum;
	}
	
	public double lossFunction(int y) {
		BinomialDist dist = new BinomialDist(y, p);
		double value = 0;
		int imin = Math.max(y - w, 0);
		for (int i = imin; i <= y; i++) {
			value += (i + w - y) * dist.prob(i);
		}
		return value;
	}
	
	public double compLossFunction(int y) {
		BinomialDist dist = new BinomialDist(y, p);
		double value = 0;
		for (int i = 0; i <= y; i++) {
			value += (y - i) * dist.prob(i);
		}
		return value;
	}
	
	public double[][] partition() {
		BinomialDist dist = new BinomialDist(ymax, p);
		double[][] result = new double[3][segNum];
		double[] p = new double[segNum];
		for (int i = 0; i < segNum; i++) {
			result[0][i] = 1.0 / (double)segNum; // 
			p[i] = result[0][i];
		}
				
		int supportUB = ymax;
		int start  = 0;
		int end = 0;
		int index = 0;
		for (int i = 0; i <= supportUB; i++) {
			if (dist.cdf(i) > p[index] + index * p[0] || i == supportUB) {
				end = i == supportUB ? i : i - 1;
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
	
	public void drawPicL() {
		double[] values = new double[ymax - ymin + 1];
		int index = 0;
		int[] support = new int[ymax - ymin + 1];
		for (int y = ymin; y < ymax; y++) {
			values[index] = lossFunction(y);
			support[index] = y;
			index++;
		}
			
		XYSeries series = new XYSeries("loss function");
		for (int k = 0; k < index; k++) {
			series.add((double) support[k], values[k]);
		}
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
	
		String title = "loss function for Binomial distribution";
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

	
	public static void main(String[] args) {
		int ymin = 30;
		int ymax = 200;
		int w = 80;
		double p = 0.1;
		int segNum = 4;
		
		BinomialPiece piece = new BinomialPiece(ymin, ymax, w, p, segNum);
//		double[][] result = piece.partition();
		piece.drawPicL();
//		System.out.println(Arrays.deepToString(result));
		

	}

}
