/**
 * 
 */
package workforce;

import java.awt.BasicStroke;
import java.awt.Color;
import java.util.Arrays;

import javax.swing.JFrame;

import org.apache.poi.xssf.usermodel.charts.AbstractXSSFChartSeries;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.openxmlformats.schemas.drawingml.x2006.main.CTRegularTextRun;

import cern.colt.function.IntIntDoubleFunction;
import milp.ComplementaryFirstOrderLossFunction;
import sdp.write.WriteToExcelTxt;
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
	int segmentNum;
	
	public TestConvexity(int w, double p, int ymax, int ymin, int segment) {
		this.w = w;
		this.p = p;
		this.ymax = ymax;
		this.ymin = ymin;
		this.segmentNum = segment;
	}
	
	public double[][] piecewise() {
		double[] slope = new double[segmentNum + 1];
		double[] intercept = new double[segmentNum + 1];
		double[] tanPointXcoe = new double[segmentNum + 1];
		double[] tanPointYcoe = new double[segmentNum + 1];
		double[] intsPointXcoe = new double[segmentNum];
		double[] intsPointYcoe = new double[segmentNum];
		double[][] result = new double[6][];
		
		int endX = ymax;
		for (int k = w; k < w*20; k ++) {
			if (Fy(k) > 0.9999) {
				endX = k;
				break;
			}
		}
		
		slope[segmentNum] = 0;
		tanPointXcoe[segmentNum] = endX;
		tanPointYcoe[segmentNum] = 0;
		intercept[segmentNum] = 0;
		for (int i = 0; i < segmentNum; i++) {
			if (i == 0) {
				slope[i] = p - 1;
				tanPointXcoe[0] = w - 1; // 切点横坐标
				tanPointYcoe[0] = (w - 1) * p + 1; // is right
				intercept[0] = w; // is right, is the y-intercept
			}
			else {
				int a = (int)tanPointXcoe[i-1];				
				for (int j = a; j <= endX; j++) {
					if (Fy(j) - Fy(a) > 1 /(double)segmentNum) {
						tanPointXcoe[i] = j;
						int b = (int)tanPointXcoe[i];
						tanPointYcoe[i] = lossFunction(b);
						slope[i] = -(1 - p)*(1 - Fy(b));
						intercept[i] = -slope[i] * tanPointXcoe[i] + tanPointYcoe[i];
						break;
					}	
				}
				if (Fy(endX) - Fy(a) <= 1 /(double)segmentNum) {
					slope[i] = slope[i-1];
					tanPointXcoe[i] = endX; // tangent point
					int b = (int)tanPointXcoe[i];
					tanPointYcoe[i] = lossFunction(b);
					intercept[i] = -slope[i] * tanPointXcoe[i] + tanPointYcoe[i];
				}
			}
		}
		for (int i = 0; i < segmentNum; i++) { // right
			intsPointXcoe[i] = tanPointYcoe[i+1] - tanPointYcoe[i] + slope[i]*tanPointXcoe[i] - slope[i+1]*tanPointXcoe[i+1];
			if (slope[i+1] > slope[i]) {
				intsPointXcoe[i] = intsPointXcoe[i] / (slope[i] -  slope[i+1]);
			}
			else {
				intsPointXcoe[i] = endX;	
			}	
			intsPointYcoe[i] = slope[i]*(intsPointXcoe[i] - tanPointXcoe[i]) + tanPointYcoe[i];
		}
		result[0] = slope;
		result[1] = intercept;
		result[2] = tanPointXcoe;
		result[3] = tanPointYcoe;
		result[4] = intsPointXcoe;
		result[5] = intsPointYcoe;
		return result;
	}
	
	/**
	 * cdf of binomial distribution
	 * @param y
	 * @return F_y(y-W)
	 */
	public double Fy(int y) {
		BinomialDist dist = new BinomialDist(y, p);
		return dist.cdf(y - w);
	}
	
	public boolean test() {
		double[] values = new double[ymax - ymin + 1];
		int index = 0;
		for (int y = ymin; y < ymax; y++) {
			values[index] = lossFunction(y);
			index++;
		}
		double[] forwardDifference = new double[ymax - ymin + 1];
		for (int i = 0; i < index; i++) {
			forwardDifference[i] = i == index ? 0 : lossFunction(i + 1 + ymin) - lossFunction(i + ymin);
		}		
		System.out.println("forward difference increasing is " + checkForwardDifference2());
		drawPicForward(forwardDifference);
		for (int i = 0; i < index-1; i++) {
			double a = forwardDifference[i];
			double b = forwardDifference[i+1] + 0.01;
			if (a > b)
				return false;
		}	
		return true;
	}
	
	/**
	 * check forward difference increasing by the transformed formula
	 * @return
	 */
	public boolean checkForwardDifference2(){
		for (int y = ymin; y < ymax; y++) {
			double value1 = 0; // g(y+1)-g(y)
			double value2 = 0; // g(y)-g(y-1)
			BinomialDist dist = new BinomialDist(y, p);
			for (int m = 0; m <= w; m++) {
				value1 +=  ((y+1)*p/(y+m-w+1)-1); // m * dist.prob(y + m - w) *
				value2 +=  (1 - (y+m-w)/(y*p)); // m * dist.prob(y + m - w) *
			}
			if (value1 < value2)
				return false;
		}
		return true;
	}
	
	public void drawPicForward(double[] forwardDifference) {
		double[] values = new double[ymax - ymin + 1];
		int index = 0;
		int[] support = new int[ymax - ymin + 1];
		for (int y = ymin; y < ymax; y++) {
			values[index] = forwardDifference[index];
			support[index] = y;
			index++;
		}
			
		XYSeries series = new XYSeries("forward difference");
		for (int k = 0; k < index; k++) {
			series.add((double) support[k], values[k]);
		}
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
		String title = "turnover rate "+  Double.toString(p) + ", " + "W = " + Integer.toString(w);
		JFreeChart chart = ChartFactory.createXYLineChart(
				title, // chart title
				"y", // x axis label
				"forward difference", // y axis label
				dataset, // data
				PlotOrientation.VERTICAL,
				true, // include legend
				false, // tooltips
				false // urls
				);
		XYPlot plot = (XYPlot)chart.getPlot();
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setSeriesLinesVisible(0, false); // 散点图，设置连线不可见
	    plot.setRenderer(renderer);
		
		ChartFrame frame = new ChartFrame("sum of binomial", chart);
		frame.pack(); // fit window to figure size
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	public void drawPicL() {
		double[] values = new double[ymax - ymin + 1];
		int index = 0;
		int[] support = new int[ymax - ymin + 1];
		double[][] gy = new double[ymax - ymin + 1][2];
		for (int y = ymin; y < ymax; y++) {
			values[index] = lossFunction(y);
			support[index] = y;
			gy[index][0] = y;
			gy[index][1] = values[index];
			index++;
		}
		
//		WriteToExcelTxt write = new WriteToExcelTxt();
//		write.writeArrayToTxt(gy, "gy.txt");
		
		double[][] result = piecewise();
		double[] slope = result[0];
		double[] intercept = result[1];
		double[] intsPointX = result[4];
		
		XYSeries series = new XYSeries("loss function");
		for (int k = 0; k < index; k++) {
			series.add((double) support[k], values[k]);
		}
		
		XYSeriesCollection dataset = new XYSeriesCollection();
		
//		dataset.addSeries(series);
		
//		XYSeries seriesLine1 = new XYSeries("line1");
//		for (double i = ymin; i <= intsPointX[0]; i = i + 0.1)
//			seriesLine1.add(i, i * slope[0] + intercept[0]);
		
		// without the below codes, it will only draw loss function
		XYSeries[] seriesLines = new XYSeries[segmentNum + 1];
		for (int k = 0; k < segmentNum + 1; k++) {
			int m = k + 1;
			seriesLines[k] = new XYSeries("line" + m);
			double b = k < segmentNum ? intsPointX[k] : ymax;
			double a = k == 0 ? ymin : intsPointX[k-1];
			for (double i = a; i <= b; i = i + 0.1)
				seriesLines[k].add(i, i * slope[k] + intercept[k]);
			dataset.addSeries(seriesLines[k]);
		}
		
		
		
		String title = "turnover rate "+  Double.toString(p) + ", " + "W = " + Integer.toString(w);
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
		renderer.setSeriesLinesVisible(0, false); // 0 is the series number
		//renderer.setSeriesPaint(1, Color.black);
		renderer.setSeriesStroke(0, new BasicStroke(50));
		renderer.setSeriesStroke(1, new BasicStroke(1));
		
		
	    plot.setRenderer(renderer);
		
		ChartFrame frame = new ChartFrame("sum of binomial", chart);
		frame.pack(); // fit window to figure size
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	/**
	 * original formulation of the loss function
	 * @param y
	 * @return
	 */
	public double lossFunction(int y) {
		BinomialDist dist = new BinomialDist(y, p);
		double value = 0;
		int imin = Math.max(y - w, 0); 
		for (int i = imin; i <= y; i++) {
			value += dist.prob(i) * (double) (i+w-y);
		}
		return value;
	}
	
	/**
	 * an equally transformation of the loss function support
	 * @param y
	 * @return
	 */
	public double lossFunctoin2(int y) {
		BinomialDist dist = new BinomialDist(y, p);
		double value = 0;
		for (int i = 0; i <= w; i++) {
			value += i * dist.prob(i+y-w);
		}
		return value;
	}
	
	public static void main(String[] args) {
		int w = 40;
		double p = 0.9;
		int ymin = w-10;
		int ymax = w*20;
		int segNum = 10;
		
		TestConvexity test = new TestConvexity(w, p, ymax, ymin, segNum);
		boolean result = test.test();
		
		System.out.println(test.lossFunction(w-1));
		
		double[][] result2 = test.piecewise();
		test.drawPicL();
		System.out.println(result);
		System.out.println(Arrays.deepToString(result2));
	}

}
