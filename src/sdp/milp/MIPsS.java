package sdp.milp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import sdp.inventory.GetPmf;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

public class MIPsS {
	double iniInventory;
	double fixedOrderingCost;
	double proportionalOrderingCost; 
	double holdingCost;
	double penaltyCost;
	int maxOrderQuantity;
	double[] demands;
	double[][][] pmf;
	boolean isForGy = false;
	
	double truncationQuantile = 0.9999;  
	double stepSize = 1; 
	double minState = -500;
	double maxState = 500;
    
	void setForGy( ) {
		isForGy = true;
	}
	
	public MIPsS(double iniInventory, double fixedOrderingCost, double proportionalOrderingCost, double holdingCost, double penaltyCost, 
			double[] demands, int maxOrderQuantity) {
		this.iniInventory = iniInventory;
		this.fixedOrderingCost = fixedOrderingCost;
		this.proportionalOrderingCost = proportionalOrderingCost;
		this.holdingCost = holdingCost;
		this.penaltyCost = penaltyCost;
		this.demands = demands;
		this.maxOrderQuantity = maxOrderQuantity;
		this.pmf = getPmf(demands);
	}
	
	double[][][] getPmf(double[] meanDemand){
		int T = demands.length;
		
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i])).toArray(PoissonDist[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
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
			return DoubleStream.iterate(0, i -> i + stepSize).limit(maxOrderQuantity + 1).toArray();
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
	  
    State stateTransition(State state, double action, boolean isForGy, double randomDemand) {
    	double nextInventory = isForGy && state.period == 1 ? state.initialInventory - randomDemand
                  : state.initialInventory + action - randomDemand;
    	nextInventory = nextInventory > maxState ? maxState : nextInventory;
    	nextInventory = nextInventory < minState ? minState : nextInventory;
    	return new State(state.period+1, nextInventory);
    }
	
    double immediateValue(State state, Double action, boolean isForGy, Double randomDemand) {
    	double fixedCost = 0, variableCost = 0, inventoryLevel = 0, holdingCosts = 0, penaltyCosts = 0;
    	if (isForGy == true && state.period == 1) {
      	  fixedCost = 0;
      	  variableCost = proportionalOrderingCost*state.initialInventory;
      	  inventoryLevel = state.initialInventory - randomDemand; 	  
        }
        else {
      	  fixedCost = action > 0 ? fixedOrderingCost : 0;
      	  variableCost = proportionalOrderingCost*action;
          inventoryLevel = state.initialInventory + action -randomDemand;	  
        }

    	holdingCosts = holdingCost*Math.max(inventoryLevel, 0);
    	penaltyCosts = penaltyCost*Math.max(-inventoryLevel, 0);
    	double totalCosts = fixedCost + variableCost + holdingCosts + penaltyCosts;
    	return totalCosts;
    }
    
    Comparator<State> keyComparator = (o1, o2) -> o1.period > o2.period ? 1 : 
		o1.period == o2.period ? o1.initialInventory > o2.initialInventory ? 1 : 
			(o1.initialInventory == o2.initialInventory ? 0 : -1) : -1;

	SortedMap<State, Double> cacheActions = new TreeMap<>(keyComparator);
	Map<State, Double> cacheValues = new HashMap<>();
	double f(State state){		
		return cacheValues.computeIfAbsent(state, s -> {
	    	  double val= Arrays.stream(s.getFeasibleActions())
	    			  .map(orderQty -> Arrays.stream(pmf[s.period-1])
	    					  .mapToDouble(p -> p[1]*immediateValue(s, orderQty, isForGy, p[0])+
	    							  (s.period < pmf.length?
	    									  p[1]*f(stateTransition(s, orderQty, isForGy, p[0])) : 0))
	    					  .sum())
	    			  .min()
	    			  .getAsDouble();
	    	  double bestOrderQty = Arrays.stream(s.getFeasibleActions())
	    			  .filter(orderQty -> Arrays.stream(pmf[s.period-1])
	    					  .mapToDouble(p -> p[1]*immediateValue(s, orderQty, isForGy, p[0])+
	    							  (s.period < pmf.length ?
	    									  p[1]*f(stateTransition(s, orderQty, isForGy, p[0])):0))
	    					  .sum() == val)
	    			  .findAny()
	    			  .getAsDouble();
	    	  cacheActions.putIfAbsent(s, bestOrderQty);
	         return val;
	      });
	}
	
	 void runSDP() {
		int period = 1;
		State initialState = new State(period, iniInventory);
		long currTime=System.currentTimeMillis();
		double finalValue = f(initialState);
		System.out.println("final optimal expected value by SDP is: " + finalValue);
		double optQ = cacheActions.get(initialState);
		System.out.println("optimal order quantity in the first priod is : " + optQ);	
		double time = (System.currentTimeMillis()-currTime)/1000;	
		System.out.println("running time is " + time + " s");
	}
	
	void drawGy(int minInventory, int maxInventory) {	
		ArrayList<Double> recordG = new ArrayList<Double>();
		XYSeries seriesG = new XYSeries("xySeries");
		setForGy();
	    for (int iniInventory = minInventory; iniInventory <= maxInventory; iniInventory = iniInventory + 1) {
	    	  int period = 1;	    	  
	    	  State initialState = new State(period, iniInventory);
	    	  double finalValue = f(initialState);
	    	  recordG.add(finalValue);
	    	  seriesG.add(iniInventory, finalValue); 	  
	    }
				
		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesG);
		
		JFreeChart chart = ChartFactory.createXYLineChart(
				"G(y) with different y", // chart title
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
	
	void drawxQ(int minInventory, int maxInventory) {	
		ArrayList<Double> recordQ = new ArrayList<Double>();
		XYSeries seriesQ = new XYSeries("xySeries");
	    for (int iniInventory = minInventory; iniInventory <= maxInventory; iniInventory = iniInventory + 1) {
	    	  int period = 1;
	    	  State initialState = new State(period, iniInventory);
	    	  f(initialState);
	    	  double optQ = cacheActions.get(initialState);
	    	  recordQ.add(optQ);
	    	  seriesQ.add(iniInventory, optQ); 	  
	    }
				
		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(seriesQ);
		
		JFreeChart chart = ChartFactory.createXYLineChart(
				"Q with different x", // chart title
				"x", // x axis label
				"Q", // y axis label
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
		
	public static void runMIPLB(double initialInventory, double[] meandemand, double fixedOrderingCost, 
			double proportionalOrderingCost, double penaltyCost, int maxOrderQuantity) {
		
	}
	
	
	
	public static void main(String[] args) {	
		double[] demands = {9, 23, 53, 29 };

		double initialInventory = 0; 
		double[] meanDemand = demands;
		double fixedOrderingCost = 500; 
		double proportionalOrderingCost = 0; 
		double holdingCost = 2;
		double penaltyCost = 10;
		int maxOrderQuantity = 100;
		
		MIPsS inventory = new MIPsS(initialInventory, fixedOrderingCost, proportionalOrderingCost, holdingCost, penaltyCost, 
				meanDemand, maxOrderQuantity);
		
		inventory.runSDP();
		int minInventory = -66;
		int maxInventory = 200;
		inventory.drawxQ(minInventory, maxInventory);
		// since comupteIfAbsent, we need initializing a new class to draw Gy; if not, java would not compute sdp again
		MIPsS inventory2 = new MIPsS(initialInventory, fixedOrderingCost, proportionalOrderingCost, holdingCost, penaltyCost, 
				meanDemand, maxOrderQuantity);
	    inventory2.drawGy(minInventory, maxInventory);
	}
}
