package sdp.inventory;

import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

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
	 * drawing a simple picture for G()
	 */
	public void drawSimpleG(double[][] yG) {
		XYSeries seriesG = new XYSeries("yQSeries");
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
				PlotOrientation.VERTICAL, false, // include legend
				true, // tooltips
				false // urls
				);

		ChartFrame frame = new ChartFrame("chen zhen'GetPmf picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	/**
	 * 
	 * drawing G() with s, S in different colors
	 */
	public void drawGAndsS(double[][] yG, double fixedOrderingCost) {
		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		XYSeries seriesG = new XYSeries("yQSeries");
		XYSeries seriesS = new XYSeries("SSeries");
		XYSeries seriesSmalls = new XYSeries("seriesSmalls");
		ArrayList<Double> recordS = new ArrayList<>();
		int N = yG.length; int SNum = 0;
		for (int i = 0; i < N; i++) {
			seriesG.add(yG[i][0], yG[i][1]);
			if (N > 2 && i > 1)
				if (yG[i - 1][1] < yG[i - 2][1] && yG[i - 1][1] < yG[i][1]) {
					seriesS.add(yG[i - 1][0], yG[i - 1][1]);
					recordS.add(yG[i - 1][1]);
					SNum++;
					// make S in descending order
					if (SNum >1) {
						if (recordS.get(SNum - 1) > recordS.get(SNum-2)) {
							seriesS.remove(SNum - 1);
							SNum--;
						}
					}
					for (int j = i - 1; j >= 0; j--) 
						if (yG[j][1] > yG[i - 1][1] + fixedOrderingCost) {
							seriesSmalls.add(yG[j][0], yG[j][1]);
							break;
						}
				}
		}

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
		// "1" is the scatter plot
		renderer.setSeriesLinesVisible(1, false);
		renderer.setSeriesShapesVisible(1, true);
		renderer.setSeriesLinesVisible(2, false);
		renderer.setSeriesShapesVisible(2, true);
		Shape circle = new Ellipse2D.Double(0, 0, 5, 5); // circle
		renderer.setSeriesShape(1, circle);
		renderer.setSeriesShape(2, circle);

		plot.setRenderer(renderer);
		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	
}
