package sdp.capacity;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.labels.StandardXYItemLabelGenerator;
import org.jfree.chart.labels.StandardXYSeriesLabelGenerator;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYItemLabelGenerator;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.util.ShapeUtilities;

import umontreal.ssj.probdist.PoissonDist;

public class CLSPforDraw {
	
	double fixedOrderingCost, proportionalOrderingCost; 
	double penaltyCost, holdingCost;
	double stepSize, minState, maxState, truncationQuantile;
	double[] demands;
	double[][][] pmf;
	int maxOrderQuantity;
	boolean isForGy = false;
     
	public CLSPforDraw(double fixedOrderingCost, double proportionalOrderingCost, double penaltyCost, 
			double holdingCost, double stepSize, double minState, double maxState, double truncationQuantile,
			double[] demands, int maxOrderQuantity) {
		this.fixedOrderingCost = fixedOrderingCost;
		this.proportionalOrderingCost = proportionalOrderingCost;
		this.holdingCost = holdingCost;
		this.penaltyCost = penaltyCost;
		this.stepSize = stepSize;
		this.minState = minState;
		this.maxState = maxState;
		this.truncationQuantile = truncationQuantile;
		this.demands = demands;
		this.maxOrderQuantity = maxOrderQuantity;
		this.pmf = getPmf(demands);
	}
	
	
	void setForGy( ) {
		isForGy = true;
	}
	
	
	double[][][] getPmf(double[] demands){
		int T = demands.length;
		PoissonDist[] distributions = IntStream.iterate(0, i -> i + 1)
				.limit(T)
				.mapToObj(i -> new PoissonDist(demands[i]))
				.toArray(PoissonDist[]::new); 
		double[] supportLB = IntStream.iterate(0, i -> i + 1)
				.limit(T)
				.mapToDouble(i -> distributions[i].inverseF(1-truncationQuantile))
				.toArray();
		double[] supportUB = IntStream.iterate(0, i -> i + 1)
				.limit(T)
				.mapToDouble(i -> distributions[i].inverseF(truncationQuantile))
				.toArray();  
		double[][][] pmf = new double[T][][];
		for (int i=0; i<T; i++)
		{
			int demandLength = (int) ((supportUB[i] - supportLB[i]+1)/stepSize);
			pmf[i] = new double[demandLength][];
			for (int j=0; j<demandLength; j++) {
				pmf[i][j] = new double[2];
				pmf[i][j][0] = supportLB[i] + j*stepSize;
				pmf[i][j][1] = distributions[i].prob((int)pmf[i][j][0])/truncationQuantile;
			}			
		}
		return pmf;
	}
	
	private class State{
		int period;
		double initialInventory;
		
		public State(int period, double initialInventory)
		{
			this.period = period;
			this.initialInventory = initialInventory;
		}
		
		public double[] getFeasibleActions(){
			return actionGenerator.apply(this);
		}
		
		@Override
		public int hashCode(){
			String hash = "";
			hash = hash + period + initialInventory;
			return hash.hashCode();
		}
		
		@Override
		public boolean equals(Object o) {
			if (o instanceof State)
				return ((State) o).period == this.period &&
						((State) o).initialInventory == this.initialInventory;
			else
				return false;
		}
		
		@Override
		public String toString() {
			return "period = " + period +", "+"initialInventory = " + initialInventory; 
		}
	}
	
	Function<State, double[]> actionGenerator = s->{
    	  return DoubleStream.iterate(0,i->i+stepSize).limit(maxOrderQuantity+1).toArray();
      };
	
    
	@FunctionalInterface
	interface StateTransitionFunction <S, A, B, R>{
		public S apply (S s, A a, B b, R r);
	}
	
	// 画 G(y) 时这个地方改变了
	StateTransitionFunction <State, Double, Boolean, Double>  stateTransition = (state, action, isForGy, randomDemand) ->{
	  	   double nextInventory = isForGy && state.period == 1 ? state.initialInventory - randomDemand
	  			                          : state.initialInventory + action - randomDemand;
	  	   nextInventory = nextInventory > maxState ? maxState : nextInventory;
	  	   nextInventory = nextInventory < minState ? minState : nextInventory;
	  	   return new State(state.period+1, nextInventory);
	     };
	
	@FunctionalInterface
	interface ImmediateValueFunction <S, A, R, B, V>{
		public V apply (S s, A a, B b, R r);
	}
	
	// 画 G(y) 时这个地方改变了
	ImmediateValueFunction<State, Double, Double, Boolean, Double> immediateValue = (state, action, isForGy, randomDemand) ->
    {
      double fixedCost = 0, variableCost = 0, inventoryLevel = 0, holdingCosts = 0, penaltyCosts = 0;
      if (isForGy == true && state.period == 1) {
    	  fixedCost = 0;
    	  variableCost = proportionalOrderingCost*state.initialInventory;
    	  inventoryLevel = state.initialInventory -randomDemand; 	  
      }
      else {
    	  fixedCost = action > 0 ? fixedOrderingCost : 0;
    	  variableCost = proportionalOrderingCost*action;
      	  inventoryLevel = state.initialInventory + action -randomDemand;	  
      }
      holdingCosts = holdingCost*Math.max(inventoryLevel, 0);
      penaltyCosts = penaltyCost*Math.max(-inventoryLevel, 0);
      double totalCosts = fixedCost + variableCost +holdingCosts + penaltyCosts;  
  	  return totalCosts;
    };
	
    
	Map<State, Double> cacheActions = new HashMap<>();
	Map<State, Double> cacheValues = new HashMap<>();
	double f(State state){
	      return cacheValues.computeIfAbsent(state, s -> {
	    	  double val= Arrays.stream(s.getFeasibleActions())
	    			  .map(orderQty -> Arrays.stream(pmf[s.period-1])
	    					  .mapToDouble(p -> p[1]*immediateValue.apply(s, orderQty, isForGy, p[0])+
	    							  (s.period < pmf.length?
	    									  p[1]*f(stateTransition.apply(s, orderQty, isForGy, p[0])) : 0))
	    					  .sum())
	    			  .min()
	    			  .getAsDouble();
	    	  double bestOrderQty = Arrays.stream(s.getFeasibleActions())
	    			  .filter(orderQty -> Arrays.stream(pmf[s.period-1])
	    					  .mapToDouble(p -> p[1]*immediateValue.apply(s, orderQty, isForGy, p[0])+
	    							  (s.period < pmf.length ?
	    									  p[1]*f(stateTransition.apply(s, orderQty, isForGy, p[0])):0))
	    					  .sum() == val)
	    			  .findAny()
	    			  .getAsDouble();
	    	  cacheActions.putIfAbsent(s, bestOrderQty);
	         return val;
	      });
	   }
	
	void drawGy(XYSeries series, XYSeries seriesS, XYSeries seriesSmalls, XYSeries seriesb) {
		//XYDataset dataset = new XYSeriesCollection(series);
		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(series);
		seriesCollection.addSeries(seriesS);
		seriesCollection.addSeries(seriesSmalls);
		seriesCollection.addSeries(seriesb);
		
		XYSeries seriesbBound = new XYSeries("seriesbBound");
		int size = seriesb.getItemCount();
		seriesbBound.add(seriesb.getX(size-1), seriesb.getY(size-1));
		seriesCollection.addSeries(seriesbBound);
		
		
		JFreeChart chart = ChartFactory.createXYLineChart(
				"Optimal cost related with initial inventory level", // chart title
				"y", // x axis label
				"G(y)", // y axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL,
				false, // include legend
				true, // tooltips
				false // urls
				);
		
		XYPlot plot = (XYPlot)chart.getPlot();
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		// "0" is the line plot
		renderer.setSeriesLinesVisible(0, true);
	    renderer.setSeriesShapesVisible(0, false);
	    
		// "1" is the scatter plot
		renderer.setSeriesLinesVisible(1, false);
		renderer.setSeriesShapesVisible(1, true);
		renderer.setSeriesLinesVisible(2, false);
		renderer.setSeriesShapesVisible(2, true);
		renderer.setSeriesLinesVisible(3, false);
		renderer.setSeriesShapesVisible(3, false);
		Shape circle =  new Ellipse2D.Double(0,0,5,5); // circle
		renderer.setSeriesShape(1, circle);
		renderer.setSeriesShape(2, circle);
		
		plot.setRenderer(renderer);
		
		NumberFormat format = NumberFormat.getNumberInstance();
		format.setMaximumFractionDigits(0); // etc.
		XYItemLabelGenerator generator1 =
		    new StandardXYItemLabelGenerator("S", format, format); //show coordinates ({1}, {2})
		XYItemLabelGenerator generator2 =
			    new StandardXYItemLabelGenerator("s", format, format); //show coordinates ({1}, {2})
		XYItemLabelGenerator generator3 =
			    new StandardXYItemLabelGenerator("({1}, {2})", format, format);
		renderer.setSeriesItemLabelGenerator(1, generator3);
		renderer.setSeriesItemLabelGenerator(2, generator2);
		//renderer.setSeriesItemLabelGenerator(3, generator3);
		renderer.setSeriesItemLabelGenerator(4, generator3);
		//renderer.setSeriesPaint(3, Color.orange);
		renderer.setBaseItemLabelsVisible(true);
		
		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	void drawxQ(XYSeries series) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);

		XYSeries orderingThreshold = new XYSeries("orderingThreshold");
		int size = series.getItemCount();
		for (int i = 0; i < size; i++) {
			if ((double)series.getY(i) > 0 && (double)series.getY(i+1) == 0)
				orderingThreshold.add(series.getX(i),series.getY(i));
		}
		dataset.addSeries(orderingThreshold);
			
		JFreeChart chart = ChartFactory.createXYLineChart(
				"Optimal Q with initial inventory level", // chart title
				"x", // x axis label
				"Q", // y axis label
				dataset, // data
				PlotOrientation.VERTICAL,
				false, // include legend
				true, // tooltips
				false // urls
				);
		XYPlot plot = (XYPlot)chart.getPlot();
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		// "0" is the line plot
		renderer.setSeriesLinesVisible(0, true);
	    renderer.setSeriesShapesVisible(0, false);
	    
		// "1" is the scatter plot
		renderer.setSeriesLinesVisible(1, false);
		renderer.setSeriesShapesVisible(1, true);
		Shape circle =  new Ellipse2D.Double(0,0,5,5); // circle
		renderer.setSeriesShape(1, circle);
		plot.setRenderer(renderer);
		
		NumberFormat format = NumberFormat.getNumberInstance();
		format.setMaximumFractionDigits(0); // etc.
		XYItemLabelGenerator generator1 =
			    new StandardXYItemLabelGenerator("({1}, {2})", format, format);
		renderer.setSeriesItemLabelGenerator(1, generator1);
		renderer.setBaseItemLabelsVisible(true);
		
		ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
	
	public static void main(String[] args) {		  
	      double[][] demands = {{20,20,20,20,20,20,20,20,20,20},
	    		  {5.4,7.2,9.6,12.2,15.4,18.6,22,25.2,28.2,30.6},
	    		  {33.2,32.4,30.6,28.2,25.2,22,18.6,15.4,12.2,9.6},
	    		  {24.2,20,15.8,14,15.8,20,24.2,26,24.2,20},
	      		  {31.4,20,8.6,4,8.6,20,31.4,36,31.4,20},
	    		  {41.8,18.2,6.6,15.8,0.4,15.2,21.8,23,44.8,4.4},
	    		  {0.4,10.2,30.4,93.4,53.6,97.8,89.2,49.6,56.2,72.6},
	    		  {9.4,16.2,47.2,78.8,32.8,57.4,101.6,78.2,150.8,138.8},
	    		  {8.8,23.2,52.8,28.8,29.2,39.6,14.8,36.6,40.8,22.8},
	    		  {9.8,37.6,12.8,55.8,90.6,44.8,44.6,103.4,58.2,109.4},
	    		  {9, 23, 53, 29, 29, 40, 15, 37, 41, 23}
	      };
	      
	      double truncationQuantile = 0.99;  
	      double stepSize = 1; 
	      double minState = -500;
	      double maxState = 300;
	      
	      // a best example for draw k convex
//	      % example: double fixedOrderingCost = 500; 
//	      % 	      double proportionalOrderingCost = 0; 
//	      % 	      double penaltyCost = 10;
//	      % 	      double[] meanDemand = {9, 23, 53, 29}; //demands[0];
//	      % 	      double holdingCost = 2;
//	      % 	      int maxOrderQuantity = 100; 
	      
	      double fixedOrderingCost = 500; 
	      double proportionalOrderingCost = 0; 
	      double penaltyCost = 10;
	      double[] meanDemand = {9, 23, 53, 29}; //demands[0];
	      double holdingCost = 2;
	      int maxOrderQuantity = 100;	//100
	      
//	      double fixedOrderingCost = 500; 
//	      double proportionalOrderingCost = 0; 
//	      double penaltyCost = 10;
//	      double[] meanDemand = {9, 23, 53, 29}; //demands[0];
//	      double holdingCost = 1;
//	      int maxOrderQuantity = 150; 
   

	      CLSPforDraw inventory = new CLSPforDraw(fixedOrderingCost, proportionalOrderingCost, penaltyCost, 
	  			holdingCost, stepSize, minState, maxState, truncationQuantile,
	  			meanDemand, maxOrderQuantity);
	      
	      ArrayList<Double> recordG = new ArrayList<Double>();
	      ArrayList<double[]> recordQ = new ArrayList<double[]>();
	      ArrayList<Double> recordS = new ArrayList<>();
	      ArrayList<Integer> recordSIndex = new ArrayList<>();
	      XYSeries seriesG = new XYSeries("xySeries");
	      XYSeries seriesQ = new XYSeries("xQSeries");
	      XYSeries seriesS = new XYSeries("SSeries");
	      XYSeries seriesSmalls = new XYSeries("seriesSmalls");
	      XYSeries seriesb = new XYSeries("seriesb");
	      int minInventory = -200;
	      int maxInventory = 200;
	      for (int initialInventory = minInventory; initialInventory <=maxInventory; initialInventory = initialInventory+1) {
	    	  int period = 1;
	    	  double[] optQ = new double[2];
	    	  State initialState = inventory.new State(period, initialInventory);
	    	  double finalValue = inventory.f(initialState);
	    	  optQ[0] = initialInventory;
	    	  optQ[1] = inventory.cacheActions.get(inventory.new State(period, initialInventory));
	    	  recordQ.add(optQ);
	    	  seriesQ.add(optQ[0], optQ[1]);
	      }
	      inventory.drawxQ(seriesQ);
	      for (int i = 0; i < recordQ.size(); i++) {
	    	  for (int j = 0; j < recordQ.get(i).length; j++) {
	    		 // System.out.printf("%.0f ", recordQ.get(i)[j]);
	    	  }
	    	  //System.out.printf("\n");
	      }
	      System.out.printf("\n\n\n\n");
	      
	      CLSPforDraw inventory2 = new CLSPforDraw(fixedOrderingCost, proportionalOrderingCost, penaltyCost, 
		  			holdingCost, stepSize, minState, maxState, truncationQuantile,
		  			meanDemand, maxOrderQuantity);
	      
	      int index = 0;
	      int Ssize = 0;
	      for (int initialInventory = minInventory; initialInventory <= maxInventory; initialInventory = initialInventory+1) {
	    	  int period = 1;
	    	  inventory2.setForGy();
	    	  State initialState = inventory2.new State(period, initialInventory);
	    	  long currTime2=System.currentTimeMillis();
	    	  double finalValue = inventory2.f(initialState);
	    	  //System.out.println("final optimal expected value is: " + finalValue);
	    	  double optQ = inventory2.cacheActions.get(inventory2.new State(period, initialInventory));
		      //System.out.println("optimal order quantity in the first priod is : " + optQ);	
	    	  double time = (System.currentTimeMillis()-currTime2)/1000;
	    	  //System.out.println("running time is " + time + " s"); 
	    	  recordG.add(finalValue);
	    	  seriesG.add(initialInventory, finalValue);  
	    	  if (recordG.size() > 2) {
	    		  if (recordG.get(index-1) < recordG.get(index-2) && recordG.get(index-1) < recordG.get(index)) {
	    			  recordS.add(recordG.get(index-1));
	    			  recordSIndex.add(index);
	    			  seriesS.add(initialInventory-1, recordG.get(index-1));
	    			  Ssize++;
	    		  }
	    	  }
	    	  if (Ssize > 1) {
	    		  if (recordS.get(Ssize - 1) > recordS.get(Ssize-2)) {
	    			  recordS.remove(Ssize - 1);
	    			  seriesS.remove(Ssize - 1);
	    			  Ssize--;
	    		  }	  
	    	  }
	    	  index++;
	      }
	      for (int i = 0; i < recordS.size(); i++) {
	    	  for (int j = recordSIndex.get(i); j >= 0; j--) {
	    		  if ((double) seriesG.getY(j) > (double) seriesS.getY(i) + fixedOrderingCost) {
	    			  seriesSmalls.add(seriesG.getX(j), seriesG.getY(j));
	    			  break;
	    		  }	  
	    	  }	    	  
//	    	  int j = recordSIndex.get(i);
//	    	  double xMinusB = (double) seriesG.getX(j)-maxOrderQuantity;
//    		  int xMinusBIndex = (int) xMinusB -minInventory + 1;
//    		  if (xMinusBIndex >=0)
//    			  seriesb.add(xMinusB, seriesG.getY(xMinusBIndex));
	      }
	      
	      
	      double lastsIndex = (double) seriesSmalls.getX(Ssize-1)-minInventory; 
	      for (int j = (int)lastsIndex; j >= 0; j--) {
	    	  double xPlusB = Math.min((double) seriesG.getX(j)+maxOrderQuantity, (double) seriesS.getX(Ssize-1));
	    	  int xPlusBIndex = (int) xPlusB -minInventory + 1;
	    	  if ((double)seriesG.getY(j) > fixedOrderingCost + (double)seriesG.getY(xPlusBIndex))
	    		  seriesb.add(seriesG.getX(j), seriesG.getY(j));  
	      }	
	      
	      inventory.drawGy(seriesG, seriesS, seriesSmalls, seriesb);	 
	      for (int i = 0; i < recordG.size(); i++) {
	    	  System.out.printf("%d %.0f\n ", minInventory + i, recordG.get(i));
	      }
	}
}
