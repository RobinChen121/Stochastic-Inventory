package workforce;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.stream.IntStream;

import javax.swing.JFrame;

import org.apache.commons.math3.analysis.function.Max;
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
	ArrayList<ArrayList<Integer>> arr = new ArrayList<>();
	int count = 0;
	
	public DemandDistribution(int[] Q, double turnoverRate) {
		this.Q = Q;
		this.turnoverRate = turnoverRate;	
	}
	
	public double[] pmf(int[] Q, double turnoverRate) {
		int sum = Arrays.stream(Q).sum();
		int T = Q.length;
		int[] support = new int[sum+1];
		double[] probs = new double[sum+1];
		BinomialDist[] dists = new BinomialDist[sum+1];
		for (int n = 1; n < sum+1; n++) {
			dists[n] = new BinomialDist(n, turnoverRate);
		}
		
		int index = 0;
		for (int n = 0; n < sum+1; n++) {
			support[n] = n;
			double prob = 0;
			int[] d = new int[T];
			
//			for (int i = 0; i < arrayList.size(); i++) {
//				d = arrayList.get(i);
//				double probj = 1;	
//				int y = 0;
//				for (int j = 0; j < T; j++) {
//					y = j == 0 ? Q[0] : y - d[j-1] + Q[j];
//					probj *= y != 0 ? dists[y].prob(d[j]) : 1;
//
//				}
//				prob += probj;
//				index++;
//			}

			
			if (T == 2) {
				for (int j = 0; j <= Math.min(n, Q[0]); j++) {
					d[0] = j; d[1] = n - j;
					double probt = 1;	
					int y = 0;
					for (int t = 0; t < T; t++) {
						y = t == 0 ? Q[0] : y - d[t-1] + Q[t];
						probt *= y != 0 ? dists[y].prob(d[t]) : 1;
					}	
					prob += probt;
//					int temp = n-j;
//					System.out.println(j + " " + temp);
					index++;
				}
			}
			else if (T == 3) {
				double p1 = 0; double p2 = 0; double p3 = 0;
				for (int i = 0; i <= Q[0]; i++) {
					d[0] = i;
					p1 = Q[0] == 0 ? 1 : dists[Q[0]].prob(i);
					int y = Q[0] - d[0] + Q[1];
					for (int j = 0; j <= y; j++) {
						d[1] = j;
						d[2] = n - i -j;
						if (n-i-j < 0)
							continue;
						p2 = y == 0 ? 1 : dists[y].prob(j);
						p3 = y - d[1] + Q[2] == 0 ? 1 : dists[y - d[1] + Q[2]].prob(n-i-j);
						prob += p1*p2*p3;
//						int temp = n-i-j;
//						System.out.println(i + " " + j + " " + temp);
						index++;
					}
				}
			}
			probs[n] = prob;
		}
		System.out.println(Arrays.stream(probs).sum());
		System.out.println(index);
		return probs;
	}
	
	/**
	 * compute pmf using Roberto's idea
	 */
	public double[] pmfR(int[] Q, double turnoverRate) {
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
		int index = 0;
		for (int k = 0; k < sum+1; k++) {
			support[k] = k;
			double prob = 0;
			int[] d = new int[T];	
		
//			for (int i = 0; i < arrayList.size(); i++) {
//				d = arrayList.get(i);
//				double probj = 1;
//				for (int j = 0; j < T; j++) {
//					probj *= dist[j].prob(d[j]);	
//				}
//				prob += probj;
//			}
			if (T == 2) {
				for (int j = 0; j <= Math.min(k, Q[0]); j++) {
					if (k-j> Q[1])
						continue;
					prob += dist[0].prob(j) * dist[1].prob(k-j);
//					int temp = k-j;
//					System.out.println(j + " " + temp);
				}
			}
			else if (T == 3) {
				for (int i = 0; i <= Math.min(k, Q[0]); i++) 
					for (int j = 0; j <= Math.min(j, Q[1]); j++) {
						if (k-i-j > Q[2] || k-i-j <0)
							continue;
						prob += dist[0].prob(j) * dist[1].prob(j)*dist[2].prob(k-i-j);
//						int temp = k-i-j;
//						System.out.println(i + " " + j + " " + temp);
						index++;
					}
			}
			
			probs[k] = prob;
		}
		System.out.println(Arrays.stream(probs).sum());
		System.out.println(index);
		return probs;
	}
	

	
	public void drawPic(double[] pmf) {
		int sum = Arrays.stream(Q).sum();
		int[] support = new int[sum+1];
		for (int k = 0; k < sum+1; k++)
			support[k] = k;
		XYSeries series = new XYSeries("pmf of the cumulative number of leaving employees");
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
	public static void main(String[] args) {
		int[] Q = {10, 10, 10};
		double turnoverRate = 0.5;
		
		DemandDistribution dist = new DemandDistribution(Q, turnoverRate);
//		int k = 3;
//		int sum = Arrays.stream(Q).sum();
//		dist.partion(2, k, k);
		
		double[] pmf = dist.pmf(Q, turnoverRate);
		dist.drawPic(pmf);
		System.out.println();
	}
	
	

}
