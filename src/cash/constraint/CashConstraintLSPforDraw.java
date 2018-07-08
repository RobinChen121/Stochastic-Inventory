package cash.constraint;

import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.DoubleStream;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.labels.StandardXYItemLabelGenerator;
import org.jfree.chart.labels.XYItemLabelGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.randvar.PoissonGen;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

/** 
* @author chen zhen 
* @version 创建时间：2018年4月15日 上午10:32:46 
* @value 类说明: 画图看 k凸性
*/
public class CashConstraintLSPforDraw {
		double price;
		double fixedOrderingCost, proportionalOrderingCost; 
		double holdingCost;
		int stepSize;
		double minCashState, maxCashState, minInventoryState, maxInventoryState;
		double truncationQuantile;
		double[] demands;
		double[][][] pmf;
		boolean isForGy = false;
		
		
		public CashConstraintLSPforDraw(double price, double fixedOrderingCost, double proportionalOrderingCost, 
				double holdingCost, int stepSize, double minInventoryState, double maxInventoryState, 
				double minCashState, double maxCashState, double truncationQuantile, double[] demands) {
			this.price = price;
			this.fixedOrderingCost = fixedOrderingCost;
			this.proportionalOrderingCost = proportionalOrderingCost;
			this.holdingCost = holdingCost;
			this.stepSize = stepSize;
			this.minInventoryState = minInventoryState;
			this.maxInventoryState = maxInventoryState;
			this.minCashState = minCashState;
			this.maxCashState = maxCashState;
			this.truncationQuantile = truncationQuantile;
			this.demands = demands;
			this.pmf = getPmf(demands);
		}
		
		void setForGy( ) {
			isForGy = true;
		}
		
		double[][][] getPmf(double[] demands){
			int T = demands.length;
			PoissonDist[] distributions = new PoissonDist[T];
			double[] supportLB = new double[T];
			double[] supportUB = new double[T];
			
			for (int i = 0; i < T; i++) {
				distributions[i] = new PoissonDist(demands[i]);
				supportLB[i] = 0; // 泊松分布特殊
				supportUB[i] = distributions[i].inverseF(truncationQuantile);
			}

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
			double iniInventory;
			double iniCash;
			
			public State(int period,  double iniCash, double iniInventory)
			{
				this.period = period;
				this.iniInventory = iniInventory;
				this.iniCash = iniCash;
			}
			
			public double[] getFeasibleActions(){
				double maxOrderingQuantity = (int) Math.max(0, (iniCash-fixedOrderingCost)/proportionalOrderingCost);
				return DoubleStream.iterate(0, i -> i + stepSize).limit((int)maxOrderingQuantity + 1).toArray();
			}
			
			@Override
			public int hashCode(){
				String hash = "";
				hash = hash + period + iniInventory + iniCash;
				return hash.hashCode();
			}
			
			@Override
			public boolean equals(Object o) {
				if (o instanceof State)
					return ((State) o).period == this.period &&
							((State) o).iniInventory == this.iniInventory &&
									((State) o).iniCash == this.iniCash;
				else
					return false;
			}
			
			@Override
			public String toString() {
				return "period = " + period +", "+"iniInventory = " + iniInventory; 
			}
		}
		
		double immediateValue(State state, double action, double randomDemand) {	
			double demand = randomDemand;
			double revenue = price*Math.min(state.iniInventory + action, demand);
			double fixedCost = action > 0 ? fixedOrderingCost : 0;
			double variableCost = proportionalOrderingCost*action;
			double inventoryLevel = state.iniInventory + action -demand;
			double holdCosts = holdingCost*Math.max(inventoryLevel, 0);
	    	
			double cashIncrement = revenue - fixedCost - variableCost - holdCosts;
	    	return cashIncrement;
	    }
		
		State stateTransition(State state, double action, double randomDemand) {
			double nextInventory = Math.max(0, state.iniInventory + action - randomDemand);
	    	nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
	    	nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
	    	
	    	double nextCash = state.iniCash + immediateValue(state, action, randomDemand);
	    	nextCash = nextCash > maxCashState ? maxCashState : nextCash;
	    	nextCash = nextCash < minCashState ? minCashState : nextCash;
	    	
	    	return new State(state.period + 1, nextCash, nextInventory);
	    }
		
		 Comparator<State> keyComparator = (o1, o2) -> o1.period > o2.period ? 1 : 
				o1.period == o2.period ? o1.iniInventory > o2.iniInventory ? 1 : 
					o1.iniInventory == o2.iniInventory ? o1.iniCash > o2.iniCash  ? 1 :
						o1.iniCash == o2.iniCash ? 0 : -1 : -1 : -1;

		SortedMap<State, Double> cacheActions = new TreeMap<>(keyComparator);
		Map<State, Double> cacheValues = new HashMap<>();
		double f(State state){		
			return cacheValues.computeIfAbsent(state, s -> {
				double[] feasibleActions = s.getFeasibleActions();
				double[][] dAndP = pmf[s.period-1]; // demandAndPossibility
				double[] QValues = new double[feasibleActions.length];
				double val = Double.MIN_VALUE;
				double bestOrderQty = 0;
				for (int i = 0; i < feasibleActions.length; i++) {
					double orderQty = feasibleActions[i];
					double thisQValue = 0;
					for (int j = 0; j < dAndP.length; j++) {
						thisQValue += dAndP[j][1]*immediateValue(s, orderQty, (int) dAndP[j][0]);
						if (s.period < pmf.length)
							thisQValue += dAndP[j][1]*f(stateTransition(s, orderQty, (int) dAndP[j][0]));
					}
					QValues[i] = thisQValue;
					if (QValues[i] > val) {
						val = QValues[i];
						bestOrderQty = orderQty;
					}
				}
				cacheActions.putIfAbsent(s, bestOrderQty);
		        return val;
		      });
		}
		
		double[][] generateSamples(double[] demands, int sampleNum){
			double[][] samples = new double[sampleNum][demands.length];
			PoissonDist[] distributions = new PoissonDist[demands.length];
			PoissonGen[] gens = new PoissonGen[demands.length];
			RandomStream genoptTable = new MRG32k3aL();
			for(int j = 0; j < demands.length; j++) {
				distributions[j] = new PoissonDist(demands[j]);
				gens[j] = new PoissonGen(genoptTable, distributions[j]);
			}
			for (int i = 0; i < sampleNum; i++) 
				for (int j = 0; j < demands.length; j++) {
					samples[i][j] = gens[j].nextInt();
			}
			return samples;
		}
		
		double simulateSamples(State iniState, Map<State, Double> cacheActions, double[][] samples) {
			double[] cashs = new double[samples.length];
			for (int i = 0; i < samples.length; i++) {
				double[] I = new double[samples[0].length];
				double[] cashBalance = new double[samples[0].length];
				double[] demand = new double[samples[0].length];
				double Q = 0, fixCost, revenue;
				State state = iniState;
				for (int t = 0; t < samples[0].length; t++)
				{
					demand[t] = samples[i][t];
					if (t==0) {
						Q = cacheActions.get(iniState); 
						I[t] = Math.max(0, iniState.iniInventory + Q - demand[t]);	
						revenue = price*Math.min(Q, demand[t]); 
					}
					else {
						try {
							Q = cacheActions.get(state); 
						}
						catch(Exception e) {
							System.out.println(t);
							}
						I[t] = Math.max(0, I[t-1] + Q - demand[t]);
						revenue = price*Math.min(Q + I[t-1], demand[t]); 
					}
					fixCost = Q > 0 ? fixedOrderingCost : 0;
					cashBalance[t] = t==0 ? iniState.iniCash + revenue - fixCost - proportionalOrderingCost*Q - holdingCost*I[t]
										  : cashBalance[t-1] + revenue - fixCost - proportionalOrderingCost*Q - holdingCost*I[t]; 
					int period = t + 1;
					state = new State(period + 1, cashBalance[t], I[t]);
				}
				cashs[i] = cashBalance[samples[0].length - 1];
			}
			return Arrays.stream(cashs).sum()/samples.length;
		}
		
		double simulatesBS(double iniCash, double[][] optsBS, double[][] samples) {
			double[] cashs = new double[samples.length];
			for (int i = 0; i < samples.length; i++) {
				double[] I = new double[samples[0].length];
				double[] cashFlow = new double[samples[0].length];
				double Q, fixCost, demand, revenue;
				for (int t = 0; t < samples[0].length; t++)
				{
					demand = samples[i][t];
					if ( t== 0) {
						Q = optsBS[t][2];
						I[t] = Math.max(0, Q - demand);
						revenue = price*Math.min(Q, demand); 
						fixCost = Q > 0 ? fixedOrderingCost : 0;
						cashFlow[t] = iniCash + revenue - fixCost- proportionalOrderingCost*Q - holdingCost*I[t];
					}
					else {
						double maxOrderQuantity = Math.max(0, (cashFlow[t-1] - fixedOrderingCost)/proportionalOrderingCost);
						if (I[t-1] < optsBS[t][0] && cashFlow[t-1] > optsBS[t][1]) 
							Q = Math.min(maxOrderQuantity, optsBS[t][2] - I[t-1]);
						else
							Q = 0;
						I[t] = Math.max(0, I[t-1] + Q - demand);
						revenue = price*Math.min(Q + I[t-1], demand); 
						fixCost = Q > 0 ? fixedOrderingCost : 0;
						cashFlow[t] = cashFlow[t-1] + revenue - fixCost- proportionalOrderingCost*Q - holdingCost*I[t];
					}
				}
				cashs[i] = cashFlow[samples[0].length - 1];
			}
			return Arrays.stream(cashs).sum()/samples.length;
		}
		
		double simulatesS(double iniCash, double[][] optsBS, double[][] samples) {
			double[] cashs = new double[samples.length];
			for (int i = 0; i < samples.length; i++) {
				double[] I = new double[samples[0].length];
				double[] cashFlow = new double[samples[0].length];
				double Q, fixCost, demand, revenue;
				for (int t = 0; t < samples[0].length; t++)
				{
					demand = samples[i][t];
					if ( t== 0) {
						Q = optsBS[t][1];
						I[t] = Math.max(0, Q - demand);
						revenue = price*Math.min(Q, demand); 
						fixCost = Q > 0 ? fixedOrderingCost : 0;
						cashFlow[t] = iniCash + revenue - fixCost- proportionalOrderingCost*Q - holdingCost*I[t];
					}
					else {
						double maxOrderQuantity = Math.max(0, (cashFlow[t-1] - fixedOrderingCost)/proportionalOrderingCost);
						if (I[t-1] >= optsBS[t][0]) 
							Q = 0;
						else
							Q = Math.min(maxOrderQuantity, optsBS[t][1] - I[t-1]);
						I[t] = Math.max(0, I[t-1] + Q - demand);
						revenue = price*Math.min(Q + I[t-1], demand); 
						fixCost = Q > 0 ? fixedOrderingCost : 0;
						cashFlow[t] = cashFlow[t-1] + revenue - fixCost- proportionalOrderingCost*Q - holdingCost*I[t];
					}
				}
				cashs[i] = cashFlow[samples[0].length - 1];
			}
			return Arrays.stream(cashs).sum()/samples.length;
		}
		
		double[][] getsBS(double iniCash, double[][] optimalTable) {
			int T = demands.length;
			double[][] optimalsS = new double[T][3];
			optimalsS[0][0] = optimalTable[0][2];
			optimalsS[0][1] = iniCash;
			optimalsS[0][2] = optimalTable[0][2] + optimalTable[0][3];
			for (int t = 1; t < T; t++) {
				final int i = t + 1;
				int recordTimes = 0;
				double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
												.toArray(double[][] :: new);
				optimalsS[t][2] = 0;
				optimalsS[t][1] = 0;
				double mark_s = 0;
				for (int j = tOptTable.length - 1; j >= 0; j--) {
					if (tOptTable[j][3] != 0) {
						if (mark_s == 0) {
							optimalsS[t][0] = tOptTable[j+1][2];
							mark_s=1;
						}
						if (tOptTable[j][2] + tOptTable[j][3] > optimalsS[t][2])
							optimalsS[t][2] = tOptTable[j][2] + tOptTable[j][3];
						recordTimes = 1;
					}
					if (tOptTable[j][3] == 0 && recordTimes == 1) {
						if (tOptTable[j][1] > optimalsS[t][1])
							optimalsS[t][1] = tOptTable[j][1];
					}		
				}
			}
			System.out.println(Arrays.deepToString(optimalsS));
			return optimalsS;
		}
		
		double[][] getsS(double iniCash, double[][] optimalTable) {
			int T = demands.length;
			double[][] optimalsS = new double[T][2];
			optimalsS[0][0] = optimalTable[0][2];
			optimalsS[0][1] = optimalTable[0][2] + optimalTable[0][3];
			for (int t = 1; t < T; t++) {
				final int i = t + 1;
				double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
												.toArray(double[][] :: new);
				for (int j = tOptTable.length - 1; j >= 0; j--) {
					if (tOptTable[j][3] != 0) {
						optimalsS[t][0] = tOptTable[j+1][2];
						optimalsS[t][1] = tOptTable[j][2] + tOptTable[j][3];
						break;
					}		
				}	
			}
			System.out.println(Arrays.deepToString(optimalsS));
			return optimalsS;
		}
		
		void drawxQ(XYSeries series) {
			XYSeriesCollection dataset = new XYSeriesCollection();
			dataset.addSeries(series);

				
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
			
			ChartFrame frame = new ChartFrame("chen zhen's picture", chart);
			frame.pack();
			frame.setVisible(true);
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		}
		
		void drawxG(XYSeries series) {
			XYSeriesCollection dataset = new XYSeriesCollection();
			dataset.addSeries(series);

				
			JFreeChart chart = ChartFactory.createXYLineChart(
					"Optimal G() with initial inventory level", // chart title
					"x", // x axis label
					"G", // y axis label
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
		
		public static void main(String[] args) {
			double[] meanDemands = {8,10,15,5};
			
			double iniCash = 50;  
			double fixedOrderingCost = 10;
			
			double proportionalOrderingCost = 1;
			double holdingCost = 5;
			double price = 10;		
			
			double truncationQuantile = 0.999;  
			int stepSize = 1;
			
			double minInventoryState = 0;
			double maxInventoryState = 200;
			double minCashState = -50;   // 会影响结果，要小于固定订货成本的负值
			double maxCashState = 1200; 
			
			CashConstraintLSPforDraw inventory = new CashConstraintLSPforDraw(price, fixedOrderingCost, proportionalOrderingCost, 
					holdingCost, stepSize, minInventoryState, maxInventoryState, 
					minCashState, maxCashState, truncationQuantile, meanDemands);
			
			double iniInventory = 0;	   
			int period = 1;
			
			State initialState = inventory.new State(period, iniCash, iniInventory);
			long currTime2=System.currentTimeMillis();
			double finalValue = inventory.f(initialState) + iniCash;
			System.out.println("final optimal expected value is: " + finalValue);
			
			double optQ = inventory.cacheActions.get(inventory.new State(period, iniCash, iniInventory));
			System.out.println("optimal order quantity in the first priod is : " + optQ);			
			double time = (System.currentTimeMillis()-currTime2)/1000;
			System.out.println("running time is " + time + " s"); 		
			
			// simulation
			Map<State, Double> cacheActions = inventory.cacheActions;
			Iterator<Map.Entry<State, Double>> iterator = cacheActions.entrySet().iterator();
			double[][] optTable = new double[cacheActions.size()][4];
			int i = 0;
			while (iterator.hasNext()) {
				Map.Entry<State, Double> entry = iterator.next();
				optTable[i++] =new double[]{entry.getKey().period, entry.getKey().iniCash, entry.getKey().iniInventory, entry.getValue()};
			}
			int sampleNum = 10000;	
			
			double[][] sBS = inventory.getsBS(iniCash, optTable);
			double simsBSFinalValue = inventory.simulatesBS(iniCash, sBS, inventory.generateSamples(meanDemands, sampleNum));
			System.out.println("final simulated sBS expected value is: " + simsBSFinalValue);
			System.out.printf("Optimality gap is: %.2f%%\n", (finalValue - simsBSFinalValue)/finalValue*100);
			
			// draw
			ArrayList<Double> recordG = new ArrayList<Double>();
		    ArrayList<double[]> recordQ = new ArrayList<double[]>();
		    XYSeries seriesQ = new XYSeries("xQSeries");
		    XYSeries seriesG = new XYSeries("xySeries");
		    int minInventory = 0; int maxInventory = 200;
		    for (int initialInventory = minInventory; initialInventory <=maxInventory; initialInventory = initialInventory+1) {
		    	  double[] optQ1 = new double[2];
		    	  State initialState2 = inventory.new State(period, iniCash, initialInventory);
				  double finalValue2 = inventory.f(initialState2); // 必须有这一步，要不 cacheActions 没有值
		    	  optQ1[0] = initialInventory;
		    	  optQ1[1] = inventory.cacheActions.get(initialState2);
		    	  recordQ.add(optQ1);
		    	  seriesQ.add(optQ1[0], optQ1[1]);
		      }
		    inventory.drawxQ(seriesQ);
		    
		    for (int initialInventory = minInventory; initialInventory <=maxInventory; initialInventory = initialInventory+1) {
		    	  inventory.setForGy();
		    	  State initialState2 = inventory.new State(period, iniCash, initialInventory);
				  double finalValue2 = inventory.f(initialState2); // 必须有这一步，要不 cacheActions 没有值
				  finalValue2 = -finalValue2;
				  recordG.add(finalValue2);
		    	  seriesG.add(initialInventory, finalValue2); 
		      }
		    inventory.drawxG(seriesG);    
		}

	}
