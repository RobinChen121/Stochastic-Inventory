package workforce;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import cern.colt.function.IntIntDoubleFunction;
import umontreal.ssj.probdist.BinomialDist;

public class DemandDistribution {
	int[] Q;
	double turnoverRate;
	ArrayList<Integer> arr = new ArrayList<>();
	
	public DemandDistribution(int[] Q, double turnoverRate) {
		this.Q = Q;
		this.turnoverRate = turnoverRate;	
	}
	
	public double[] pmf(int[] Q, double turnoverRate) {
		int sum = Arrays.stream(Q).sum();
		int[] support = new int[sum+1];
		double[] probs = new double[sum+1];
		int T = Q.length;
		double[] p = new double[T];
		BinomialDist[] dist = new BinomialDist[T];
		for (int i = 0; i < T; i++) {
			p[i] = 1 - Math.pow(1-turnoverRate, T - i);
			dist[i] = new BinomialDist(Q[i], p[i]);
		}
		for (int k = 0; k < sum+1; k++) {
			support[k] = k;
			double prob = 0;
			double prob1 = 0;
			double prob2 = 0;
			for (int i = 0; i < k; i++) {
				prob1 = dist[0].prob(i);
				prob2 = dist[1].prob(k-i);
				prob += prob1*prob2;
			}
			probs[k] = prob;
		}
		return probs;
	}
	
	public void drawPic(int[] support, double[] pmf) {
		int sum = Arrays.stream(Q).sum();
		XYSeries series = new XYSeries("xySeries");
		for (int i = 0; i < sum + 1; i++) {
			series.add((double) support[i], pmf[i]);
		}
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
		String title = "turnover rate "+  Double.toString(turnoverRate) + ", Q = " + Arrays.toString(Q);
		JFreeChart chart = ChartFactory.createXYLineChart(
				title, // chart title
				"k", // x axis label
				"pmf", // y axis label
				dataset, // data
				PlotOrientation.VERTICAL,
				false, // include legend
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
	
	public void partion(int n, int k) {		
		if (k == 0)
			return;
		if (k == 1) {
			System.out.println(Integer.toString(n));
		}
		
		for (int i = 0; i <= n; i++) {
			System.out.print(Integer.toString(i));
			System.out.print(" ");
			int newN = n-i;
			int newK = k-1;
			partion(n-i, k-1);
		}
	}
	
	public static void main(String[] args) {
		int[] Q = {50, 40};
		double turnoverRate = 0.6;
		
		DemandDistribution dist = new DemandDistribution(Q, turnoverRate);
		double[] pmf = dist.pmf(Q, turnoverRate);
		dist.partion(3, 3);

		
	}
	
	

}
