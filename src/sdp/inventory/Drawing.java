package sdp.inventory;

import java.awt.Color;
import java.awt.Font;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.text.NumberFormat;
import java.util.ArrayList;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.labels.StandardXYItemLabelGenerator;
import org.jfree.chart.labels.XYItemLabelGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.RectangleEdge;


/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 10, 2018---6:54:43 PM
 * @description: drawing pictures For SDP problems to analysis
 */

public class Drawing {

	/**
	 * drawing a picture about optimal ordering quantities for different initial
	 * inventory levels
	 */
	public void drawXQ(double[][] xQ) {
		XYSeries seriesQ = new XYSeries("xQSeries");
		int N = xQ.length;
		for (int i = 0; i < N; i++) {
			seriesQ.add(xQ[i][0], xQ[i][1]);
		}

		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesQ);

		JFreeChart chart = ChartFactory.createXYLineChart("Optimal Q with different x", // chart title
				"x", // x axis label
				"Q", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, false, // include legend
				true, // tooltips
				false // urls
				);
		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	/**
	 * 
	 * drawing a simple picture for G(), with initial cash
	 */
	public void drawSimpleG(double[][] yG, double iniCash) {
		XYSeries seriesG = new XYSeries("yQSeries");
		int N = yG.length;
		for (int i = 0; i < N; i++) {
			seriesG.add(yG[i][0], yG[i][1]);
		}

		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesG);
		String cashString = String.valueOf(iniCash);
		String title = "G(y) with different order-up-to level y, B0 = " + cashString;

		JFreeChart chart = ChartFactory.createXYLineChart(title, // chart title
				"y", // x axis label
				"G(y)", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, false, // include legend
				true, // tooltips
				false // urls
				);

		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	/**
	 * 
	 * drawing a simple picture for G(), with initial cash and a String input
	 */
	public void drawSimpleG(double[][] yG, double iniCash, String str) {
		XYSeries seriesG = new XYSeries("yQSeries");
		int N = yG.length;
		for (int i = 0; i < N; i++) {
			seriesG.add(yG[i][0], yG[i][1]);
		}

		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesG);
		String cashString = String.valueOf(iniCash);
		String title = str + "(y) with different y, B0 = " + cashString;

		JFreeChart chart = ChartFactory.createXYLineChart(title, // chart title
				"y", // x axis label
				"G(y)", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, false, // include legend
				true, // tooltips
				false // urls
				);

		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	/**
	 * 
	 * drawing G() with s, S in different colors
	 */
	public static void drawGAndsS(double[][] yG, double fixedOrderingCost) {
		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		XYSeries seriesG = new XYSeries("yQSeries");
		XYSeries seriesS = new XYSeries("SSeries");
		XYSeries seriesSmalls = new XYSeries("seriesSmalls");
		ArrayList<Double> recordS = new ArrayList<>();
		ArrayList<Double> records = new ArrayList<>();
		int N = yG.length; int SNum = 0;
		for (int i = 0; i < N; i++) {
			seriesG.add(yG[i][0], yG[i][1]);
		}
		
		for (int i = 0; i <= N - 1; i++) {	
			if (N > 2 && i > 1) {
				if ( yG[i - 1][1] < yG[i - 2][1] - 0.01 && yG[i - 1][1] < yG[i][1] - 0.01) {
					seriesS.add(yG[i - 1][0], yG[i - 1][1]);
					recordS.add(yG[i - 1][1]);
					SNum++;
					// make S in descending order
					if (SNum >1) {
						if (recordS.get(SNum - 1) > recordS.get(SNum-2)) {
							seriesS.remove(SNum - 1);
							recordS.remove(SNum - 1);
							//seriesSmalls.remove(SNum - 1);
							//records.remove(SNum - 1);
							SNum--;
						}
					}
					for (int j = i - 1; j >= 0; j--) 
						if (yG[j][1] > yG[i - 1][1] + fixedOrderingCost) {
							seriesSmalls.add(yG[j][0], yG[j][1]);
							records.add(yG[i - 1][1]);
							System.out.printf("the slope at s is: %.2f\n", yG[j][1] - yG[j - 1][1]);
							break;
						}
				}
				if (i == N - 1 && SNum == 0) {
					seriesS.add(yG[i][0], yG[i][1]);
					for (int j = i - 1; j >= 0; j--) 
						if (yG[j][1] > yG[i - 1][1] + fixedOrderingCost) {
							seriesSmalls.add(yG[j][0], yG[j][1]);
							break;
						}
				}
			}
		}
		
		// show the slope below s
		seriesCollection.addSeries(seriesG);
		seriesCollection.addSeries(seriesS);
		seriesCollection.addSeries(seriesSmalls);
		JFreeChart chart = ChartFactory.createXYLineChart("Optimal cost related with initial inventory level", // chart
				// title
				"y", // x axis label
				"G(y)", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, false, // include legend
				true, // tooltips
				false // urls
				);

		XYPlot plot = (XYPlot) chart.getPlot();
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();

		// "0" is the line plot
		renderer.setSeriesLinesVisible(0, true);
		renderer.setSeriesShapesVisible(0, false);
		// "1" and "2" is the scatter plot
		renderer.setSeriesLinesVisible(1, false);
		renderer.setSeriesShapesVisible(1, true);
		renderer.setSeriesLinesVisible(2, false);
		renderer.setSeriesShapesVisible(2, true);
		Shape circle = new Ellipse2D.Double(0, 0, 5, 5); // circle
		renderer.setSeriesShape(1, circle);
		renderer.setSeriesShape(2, circle);
		
		NumberFormat format = NumberFormat.getNumberInstance();
		format.setMaximumFractionDigits(0); 
		XYItemLabelGenerator generator1 = new StandardXYItemLabelGenerator("S ({1}, {2})", format, format); 
		XYItemLabelGenerator generator2 = new StandardXYItemLabelGenerator("s ({1}, {2})", format, format); 
		//XYItemLabelGenerator generator3 = new StandardXYItemLabelGenerator("({1}, {2})", format, format); // coordinates
		renderer.setSeriesItemLabelGenerator(1, generator1);
		renderer.setSeriesItemLabelGenerator(2, generator2);
		
		renderer.setBaseItemLabelsVisible(true);

		plot.setRenderer(renderer);
		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	/**
	 * 
	 * drawing a picture for C() and different B with fixed x
	 */
	public void drawBC(double[][] BC) {
		XYSeries seriesG = new XYSeries("BCSeries");
		int N = BC.length;
		for (int i = 0; i < N; i++) {
			seriesG.add(BC[i][0], BC[i][1]);
		}

		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesG);

		JFreeChart chart = ChartFactory.createXYLineChart("C() with different ini cash B", // chart title
				"B", // x axis label
				"C()", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, false, // include legend
				true, // tooltips
				false // urls
				);

		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	/**
	 * 
	 * drawing a picture for C() and different x
	 */
	public static void drawXC(double[][] XC) {
		XYSeries seriesG = new XYSeries("BCSeries");
		int N = XC.length;
		for (int i = 0; i < N; i++) {
			seriesG.add(XC[i][0], XC[i][1]);
		}

		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesG);

		JFreeChart chart = ChartFactory.createXYLineChart("C() with different ini cash X", // chart title
				"X", // x axis label
				"C()", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, false, // include legend
				true, // tooltips
				false // urls
				);

		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	
	/**
	 * 
	 * drawing a picture for Q() and different B with fixed x
	 */
	public void drawBQ(double[][] BQ) {
		XYSeries seriesG = new XYSeries("BQSeries");
		int N = BQ.length;
		for (int i = 0; i < N; i++) {
			seriesG.add(BQ[i][0], BQ[i][1]);
		}

		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesG);

		JFreeChart chart = ChartFactory.createXYLineChart("Q with different ini cash B", // chart title
				"B", // x axis label
				"Q", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, false, // include legend
				true, // tooltips
				false // urls
				);

		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	
	/**
	 * 
	 * drawing a simple picture for G()
	 */	
	public static void drawSimpleG(double[][] yG) {
		XYSeries seriesG = new XYSeries("yGSeries");
		int N = yG.length;
		for (int i = 0; i < N; i++) {
			seriesG.add(yG[i][0], yG[i][1]);
		}

		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesG);

		JFreeChart chart = ChartFactory.createXYLineChart("G(y) with different order-up-to level y", // chart title
				"y", // x axis label
				"G(y)", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, 
				false, // include legend
				true, // tooltips
				false // urls
				);

		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	/**
	 * 
	 * drawing a picture for GA and GB for fixed initial cash R
	 */
	public void drawTwoG(double[][] GA, double[][] GB, double iniCash) {
		XYSeries seriesGA = new XYSeries("GA");
		XYSeries seriesGB = new XYSeries("GB");
		int N = GA.length;
		for (int i = 0; i < N; i++) {
			seriesGA.add(GA[i][0], GA[i][1]);
			seriesGB.add(GB[i][0], GB[i][1]);
		}

		XYSeriesCollection seriesCollectionA = new XYSeriesCollection(seriesGA);
		XYSeriesCollection seriesCollectionB = new XYSeriesCollection(seriesGB);
		
		String cashString = String.valueOf(iniCash);
		String title = "G(y) with different order-up-to level y, R0 = " + cashString;
		
		JFreeChart chart = ChartFactory.createXYLineChart(title, // chart title
				"y", // x axis label
				"G()", // y axis label
				seriesCollectionA, // data
				PlotOrientation.VERTICAL, 
				false, //true, // include legend
				true, // tooltips
				false // urls
				);
		
//		LegendTitle legend = chart.getLegend();
//		legend.setPosition(RectangleEdge.TOP);		
		
		XYPlot plot = (XYPlot) chart.getPlot();
		
		// "0" is the GA plot, "1" is the GB plot		
		plot.setDataset(0, seriesCollectionA); // first data set
		plot.setDataset(1, seriesCollectionB); // second data set
//		//plot.mapDatasetToRangeAxis(1, 0); // same axis, different data set
		
		XYItemRenderer renderer0 = new XYLineAndShapeRenderer(true, false); // boolean lines, boolean shapes
		XYItemRenderer renderer1 = new XYLineAndShapeRenderer(true, false); // boolean lines, boolean shapes
		plot.setRenderer(0, renderer0); 
		plot.setRenderer(1, renderer1); 
//		
		//plot.getRendererForDataset(plot.getDataset(0)).setSeriesPaint(0, null); 
		//plot.getRendererForDataset(plot.getDataset(1)).setSeriesPaint(1, null);
		
		// legend
		LegendTitle lt = new LegendTitle(plot);
		lt.setItemFont(new Font("Dialog", Font.PLAIN, 15));
		lt.setBackgroundPaint(new Color(200, 200, 255, 100));
		lt.setFrame(new BlockBorder(Color.white));
		lt.setPosition(RectangleEdge.BOTTOM);
		XYTitleAnnotation ta = new XYTitleAnnotation(0.98, 0.95, lt,RectangleAnchor.TOP_RIGHT);

		ta.setMaxWidth(1.5);
		plot.addAnnotation(ta);
		
		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	/**
	 * 
	 * drawing a picture for GA and GB for fixed initial inventory and different R
	 */
	public void drawTwoGR(double[][] GA, double[][] GB, double iniInventory) {
		XYSeries seriesGA = new XYSeries("GA");
		XYSeries seriesGB = new XYSeries("GB");
		int N = GA.length;
		for (int i = 0; i < N; i++) {
			seriesGA.add(GA[i][0], GA[i][1]);
			seriesGB.add(GB[i][0], GB[i][1]);
		}

		XYSeriesCollection seriesCollectionA = new XYSeriesCollection(seriesGA);
		XYSeriesCollection seriesCollectionB = new XYSeriesCollection(seriesGB);
		
		String cashString = String.valueOf(iniInventory);
		String title = "G(y) with different initial cash R, y0 = " + cashString;
		
		JFreeChart chart = ChartFactory.createXYLineChart(title, // chart title
				"R", // x axis label
				"G()", // y axis label
				seriesCollectionA, // data
				PlotOrientation.VERTICAL, 
				false, //true, // include legend
				true, // tooltips
				false // urls
				);
			
		
		XYPlot plot = (XYPlot) chart.getPlot();
		
		// "0" is the GA plot, "1" is the GB plot		
		plot.setDataset(0, seriesCollectionA); // first data set
		plot.setDataset(1, seriesCollectionB); // second data set
		
		XYItemRenderer renderer0 = new XYLineAndShapeRenderer(true, false); // boolean lines, boolean shapes
		XYItemRenderer renderer1 = new XYLineAndShapeRenderer(true, false); // boolean lines, boolean shapes
		plot.setRenderer(0, renderer0); 
		plot.setRenderer(1, renderer1); 
		
		// legend
		LegendTitle lt = new LegendTitle(plot);
		lt.setItemFont(new Font("Dialog", Font.PLAIN, 15));
		lt.setBackgroundPaint(new Color(200, 200, 255, 100));
		lt.setFrame(new BlockBorder(Color.white));
		lt.setPosition(RectangleEdge.BOTTOM);
		XYTitleAnnotation ta = new XYTitleAnnotation(0.98, 0.95, lt,RectangleAnchor.TOP_RIGHT);

		ta.setMaxWidth(1.5);
		plot.addAnnotation(ta);
		
		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	
	/**
	 * @param GA
	 * @param GB
	 * @param iniCash
	 * @return the intersection point of GA and GB
	 * @date: May 18, 2020, 10:38:49 AM 
	 */
	public double[] intersectionPoint(double[][] GA, double[][] GB, double iniCash) {
		int L = GA.length;
		double[] point = new double[3];
		for (int i = 0; i < L; i++) {
			if (GA[i][1] - GB[i][1] < 0.1) {
				point[0] = GA[i][0];
				point[1] = iniCash;
				point[2] = GA[i][1];
				break;
			}
		}
		return point;
	}	
}
