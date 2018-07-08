package variable.capacity;

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
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import umontreal.ssj.probdist.PoissonDist;

public class CLSPforDraw {
	
	double fixedOrderingCost, proportionalOrderingCost; 
	double penaltyCost, holdingCost;
	double stepSize, minState, maxState, truncationQuantile;
	double[] demands;
	double[][][] pmf;
	int maxOrderQuantity;
     
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
	
	void draw(XYSeries series) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
		JFreeChart chart = ChartFactory.createXYLineChart(
				"Optimal cost related with initial inventory level", // chart title
				"y", // x axis label
				"G(y)", // y axis label
				dataset, // data
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
	
	void drawS(XYSeries series) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
		JFreeChart chart = ChartFactory.createXYLineChart(
				"Optimal cost related with initial inventory level", // chart title
				"y", // x axis label
				"G(y)", // y axis label
				dataset, // data
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
	interface StateTransitionFunction <S, A, R>{
		public S apply (S s, A a, R r);
	}
	StateTransitionFunction <State, Double, Double>  stateTransition = (state, action, randomDemand) ->{
	  	   double nextInventory = state.initialInventory + action - randomDemand;
	  	   nextInventory = nextInventory > maxState ? maxState : nextInventory;
	  	   nextInventory = nextInventory < minState ? minState : nextInventory;
	  	   return new State(state.period+1, state.initialInventory + action - randomDemand);
	     };
	
	@FunctionalInterface
	interface ImmediateValueFunction <S, A, R, V>{
		public V apply (S s, A a, R r);
	}
	ImmediateValueFunction<State, Double, Double, Double> immediateValue = (state, action, randomDemand) ->
    {
  	  double fixedCost = action > 0 ? fixedOrderingCost : 0;
  	  double variableCost = proportionalOrderingCost*action;
  	  double inventoryLevel = state.initialInventory + action -randomDemand;
  	  double holdingCosts = holdingCost*Math.max(inventoryLevel, 0);
  	  double penaltyCosts = penaltyCost*Math.max(-inventoryLevel, 0);
  	  double totalCosts = fixedCost + variableCost +holdingCosts + penaltyCosts;
  	  return totalCosts;
    };
	
    
	Map<State, Double> cacheActions = new HashMap<>();
	Map<State, Double> cacheValues = new HashMap<>();
	double f(State state){
	      return cacheValues.computeIfAbsent(state, s -> {
	    	  double val= Arrays.stream(s.getFeasibleActions())
	    			  .map(orderQty -> Arrays.stream(pmf[s.period-1])
	    					  .mapToDouble(p -> p[1]*immediateValue.apply(s, orderQty, p[0])+
	    							  (s.period < pmf.length?
	    									  p[1]*f(stateTransition.apply(s, orderQty, p[0])) : 0))
	    					  .sum())
	    			  .min()
	    			  .getAsDouble();
	    	  double bestOrderQty = Arrays.stream(s.getFeasibleActions())
	    			  .filter(orderQty -> Arrays.stream(pmf[s.period-1])
	    					  .mapToDouble(p -> p[1]*immediateValue.apply(s, orderQty, p[0])+
	    							  (s.period < pmf.length ?
	    									  p[1]*f(stateTransition.apply(s, orderQty, p[0])):0))
	    					  .sum() == val)
	    			  .findAny()
	    			  .getAsDouble();
	    	  cacheActions.putIfAbsent(s, bestOrderQty);
	         return val;
	      });
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
	      double minState = -1000;
	      double maxState = 500;
	      
	      double fixedOrderingCost = 500; 
	      double proportionalOrderingCost = 0; 
	      double penaltyCost = 10;
	      double[] meanDemand = {9, 23, 53, 29}; //demands[0];
	      double holdingCost = 2;
	      int maxOrderQuantity = 60;	      

	      CLSPforDraw inventory = new CLSPforDraw(fixedOrderingCost, proportionalOrderingCost, penaltyCost, 
	  			holdingCost, stepSize, minState, maxState, truncationQuantile,
	  			meanDemand, maxOrderQuantity);
	      
	      ArrayList<Double> recordGy = new ArrayList<Double>();
	      ArrayList<Double> recordQ = new ArrayList<Double>();
	      XYSeries series = new XYSeries("xySeries");
	      XYSeries series2 = new XYSeries("xQSeries");
	      for (int initialInventory = -200; initialInventory <=200; initialInventory = initialInventory+1) {
	    	  int period = 1;
	    	  State initialState = inventory.new State(period, initialInventory);
	    	  long currTime2=System.currentTimeMillis();
	    	  double finalValue = inventory.f(initialState);
	    	  System.out.println("final optimal expected value is: " + finalValue);
	    	  double optQ = inventory.cacheActions.get(inventory.new State(period, initialInventory));
		      System.out.println("optimal order quantity in the first priod is : " + optQ);	
	    	  double time = (System.currentTimeMillis()-currTime2)/1000;
	    	  System.out.println("running time is " + time + " s"); 
	    	  recordGy.add(finalValue);
	    	  recordQ.add(optQ);
	    	  series.add(initialInventory, finalValue);
	    	  series2.add(initialInventory, optQ);
	      }
	      inventory.drawS(series);
	      inventory.draw(series2);
	}
}
