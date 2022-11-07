/**
 * 
 */
package workforce;

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

import milp.ComplementaryFirstOrderLossFunction;
import umontreal.ssj.probdist.BinomialDist;

/**
 * @author chen
 * @date 2022 Oct 27, 15:30 pm
 * @describe: test the convexity of loss function
 *
 */
public class TestConvexity {
	int w;
	double p;
	int ymax;
	int ymin;
	
	public TestConvexity(int w, double p, int ymax, int ymin) {
		this.w = w;
		this.p = p;
		this.ymax = ymax;
		this.ymin = ymin;
	}

	public boolean test() {
		double[] values = new double[ymax - ymin + 1];
		int index = 0;
		for (int y = ymin; y < ymax; y++) {
			values[index] = lossFunctoin(y);
			index++;
		}
		double[] forwardDifference = new double[ymax - ymin + 1];
		for (int i = 0; i < index; i++) {
			forwardDifference[i] =  i == index ? 0 : lossFunctoin(i + 1 + ymin) - lossFunctoin(i + ymin);
		}
		for (int i = 0; i < index-1; i++) {
			double a = forwardDifference[i];
			double b = forwardDifference[i+1] + 0.01;
			if (a > b)
				return false;
		}
		return true;
	}
	
	public void drawPic() {
		double[] values = new double[ymax - ymin + 1];
		int index = 0;
		int[] support = new int[ymax - ymin + 1];
		for (int y = ymin; y < ymax; y++) {
			values[index] = lossFunctoin(y);
			support[index] = y;
			index++;
		}
			
		XYSeries series = new XYSeries("loss function");
		for (int k = 0; k < index; k++) {
			series.add((double) support[k], values[k]);
		}
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
		String title = "turnover rate "+  Double.toString(p) + ", " + "W = " + Integer.toString(w);
		JFreeChart chart = ChartFactory.createXYLineChart(
				title, // chart title
				"k", // x axis label
				"pmf", // y axis label
				dataset, // data
				PlotOrientation.VERTICAL,
				true, // include legend
				false, // tooltips
				false // urls
				);
		XYPlot plot = (XYPlot)chart.getPlot();
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesLinesVisible(0, true); // 设置连线不可见
	        plot.setRenderer(renderer);
		
		ChartFrame frame = new ChartFrame("sum of binomial", chart);
		frame.pack(); // fit window to figure size
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	public double lossFunctoin(int y) {
		BinomialDist dist = new BinomialDist(y, p);
		double value = 0;
		for (int i = Math.max(y - w, 0); i <= y; i++) {
			value += dist.prob(i) * (double) (i + w - y);
		}
		return value;
	}
	
	public double[] lossFunctoin2(int y) {
		BinomialDist dist = new BinomialDist(y, p);
		double[] values = new double[y + 1];
		for (int i = 0; i <= y; i++) {
			values[i] = i * dist.prob(i);
		}
		return values;
	}
	
	public static void main(String[] args) {
		int w = 40;
		double p = 0.9;
		int ymin = w - 30;
		int ymax = w + 100;
		
		TestConvexity testConvexity = new TestConvexity(w, p, ymax, ymin);
		boolean result = testConvexity.test();
		
		testConvexity.drawPic();
		double[] values = testConvexity.lossFunctoin2(100);
		System.out.println(result);

	}

}
