package variable.capacity;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.DoubleStream;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import sdp.sampling.SampleFactory;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.randvar.PoissonGen;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

public class StructureDPforThreesS {
	double fixedOrderingCost, proportionalOrderingCost; 
	double penaltyCost, holdingCost;
	double stepSize, minState, maxState, truncationQuantile;
	double[] demands;
	double[][][] pmf;
	int[] maxOrderQuantity;
     
	public StructureDPforThreesS(double fixedOrderingCost, double proportionalOrderingCost, double penaltyCost, 
			double holdingCost, double stepSize, double minState, 
			double maxState, double truncationQuantile,
			double[] demands, int[] maxOrderQuantity) {
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
		PoissonDist[] distributions = new PoissonDist[T];
		double[] supportLB = new double[T];
		double[] supportUB = new double[T];
		
		for (int i = 0; i < T; i++) {
			distributions[i] = new PoissonDist(demands[i]);
			supportLB[i] = 0;//distributions[i].inverseF(1-truncationQuantile);
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
		double initialInventory;
		
		public State(int period, double initialInventory)
		{
			this.period = period;
			this.initialInventory = initialInventory;
		}
		
		public double[] getFeasibleActions(){
			return DoubleStream.iterate(0, i -> i + stepSize).limit(maxOrderQuantity[period-1] + 1).toArray();
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
	  
    State stateTransition(State state, double action, double randomDemand) {
    	double nextInventory = state.initialInventory + action - randomDemand;
    	nextInventory = nextInventory > maxState ? maxState : nextInventory;
    	nextInventory = nextInventory < minState ? minState : nextInventory;
    	return new State(state.period+1, state.initialInventory + action - randomDemand);
    }
	
    double immediateValue(State state, Double action, Double randomDemand) {
    	double fixedCost = action > 0 ? fixedOrderingCost : 0;
    	double variableCost = proportionalOrderingCost*action;
    	double inventoryLevel = state.initialInventory + action -randomDemand;
    	double holdingCosts = holdingCost*Math.max(inventoryLevel, 0);
    	double penaltyCosts = penaltyCost*Math.max(-inventoryLevel, 0);
    	double totalCosts = fixedCost + variableCost +holdingCosts + penaltyCosts;
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
	    					  .mapToDouble(p -> p[1]*immediateValue(s, orderQty, p[0])+
	    							  (s.period < pmf.length?
	    									  p[1]*f(stateTransition(s, orderQty, p[0])) : 0))
	    					  .sum())
	    			  .min()
	    			  .getAsDouble();
	    	  double bestOrderQty = Arrays.stream(s.getFeasibleActions())
	    			  .filter(orderQty -> Arrays.stream(pmf[s.period-1])
	    					  .mapToDouble(p -> p[1]*immediateValue(s, orderQty, p[0])+
	    							  (s.period < pmf.length ?
	    									  p[1]*f(stateTransition(s, orderQty, p[0])):0))
	    					  .sum() == val)
	    			  .findAny()
	    			  .getAsDouble();
	    	  cacheActions.putIfAbsent(s, bestOrderQty);
	         return val;
	      });
	}
	
	double[][] generateSamples(double[] demands, int sampleNum){
		double[][] samples = new double[sampleNum][demands.length];
		PoissonDist[] distributions = new PoissonDist[demands.length];
		PoissonGen[] gens = new PoissonGen[demands.length];
		RandomStream genArr = new MRG32k3aL();
		for(int j = 0; j < demands.length; j++) {
			distributions[j] = new PoissonDist(demands[j]);
			gens[j] = new PoissonGen(genArr, distributions[j]);
		}
		for (int i = 0; i < sampleNum; i++) 
			for (int j = 0; j < demands.length; j++) {
				samples[i][j] = gens[j].nextInt();
		}
		return samples;
	}
	
	double simulateSamples(double[][] optsS, double[][] samples) {
		double[] costs = new double[samples.length];
		for (int i = 0; i < samples.length; i++) {
			double[] I = new double[samples[0].length];
			double Q, fixCost, Iplus, Iminus, demand;
			for (int t = 0; t < samples[0].length; t++)
			{
				demand = samples[i][t];
				if ( t== 0) {
					Q = optsS[t][1];
					I[t] = Q - demand;
				}
				else {
					if (I[t-1] < optsS[t][0])
						Q = Math.min(maxOrderQuantity[t], optsS[t][1] - I[t-1]);
					else if (optsS[t][0] <= I[t-1] && I[t-1] < optsS[t][2])
						Q = Math.min(maxOrderQuantity[t], optsS[t][3] - I[t-1]);
					else if (optsS[t][2] <= I[t-1] && I[t-1] < optsS[t][4])
						Q = Math.min(maxOrderQuantity[t], optsS[t][5] - I[t-1]);
					else
						Q = 0.0;
					I[t] = I[t-1] + Q - demand;
				}
				Iplus = I[t] >= 0 ? I[t] : 0;
				Iminus = I[t] < 0 ? -I[t] : 0;
				fixCost = Q > 0 ? fixedOrderingCost : 0;
				costs[i] += fixCost + proportionalOrderingCost * Q + holdingCost * Iplus + penaltyCost * Iminus;
			}
		}
		return Arrays.stream(costs).sum()/samples.length;
	}
	
	int[] levelNum(double[][] optTable) {
		ArrayList<Integer> indexArr = new ArrayList<>();
		boolean mark = false;
		for (int j = 0; j < optTable.length; j++) {
			if (optTable[j][2] < maxOrderQuantity[(int) optTable[j][0]-1] && !mark) {
				mark = true;
			}
			else if (optTable[j][2] == maxOrderQuantity[(int) optTable[j][0]-1] && mark) {
				mark = false;
				indexArr.add(j);
			}
			if (optTable[j][2] == 0) {
				indexArr.add(j);
				break;
			}
		}
		return indexArr.stream().mapToInt(p -> p.intValue()).toArray();
	}
    
	double[][] getsS(double[][] optimalTable) {
		int T = demands.length;
		double[][] optimalsS = new double[T][6];
		optimalsS[0][0] = optimalTable[0][1];
		optimalsS[0][1] = optimalTable[0][1] + optimalTable[0][2];
		optimalsS[0][2] = optimalTable[0][1];
		optimalsS[0][3] = optimalTable[0][1] + optimalTable[0][2];
		optimalsS[0][4] = optimalTable[0][1];
		optimalsS[0][5] = optimalTable[0][1] + optimalTable[0][2];
		for (int t = 1; t < T; t++) {
			final int i = t + 1;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
											.toArray(double[][] :: new);
			int[] numIndex = levelNum(tOptTable);
			if (numIndex.length == 1 && numIndex[0] != 0) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				optimalsS[t][2] = tOptTable[numIndex[0]][1];
				optimalsS[t][3] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				optimalsS[t][4] = tOptTable[numIndex[0]][1];
				optimalsS[t][5] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
			}
			else if (numIndex.length == 1 && numIndex[0] == 0) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0]][1];
				optimalsS[t][2] = tOptTable[numIndex[0]][1];
				optimalsS[t][3] = tOptTable[numIndex[0]][1];
				optimalsS[t][4] = tOptTable[numIndex[0]][1];
				optimalsS[t][5] = tOptTable[numIndex[0]][1];
			}
			else if (numIndex.length == 0) {
				optimalsS[t][0] = tOptTable[tOptTable.length - 1][1];
				optimalsS[t][1] = maxOrderQuantity[t-1] * 10;
				optimalsS[t][2] = tOptTable[tOptTable.length - 1][1];
				optimalsS[t][3] = maxOrderQuantity[t-1] * 10;
			}	
			else if (numIndex.length == 2) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				optimalsS[t][2] = tOptTable[numIndex[1]][1];
				optimalsS[t][3] = tOptTable[numIndex[1] - 1][1] + tOptTable[numIndex[1] - 1][2];
				optimalsS[t][4] = tOptTable[numIndex[1]][1];
				optimalsS[t][5] = tOptTable[numIndex[1] - 1][1] + tOptTable[numIndex[1] - 1][2];
			}
			else if (numIndex.length == 3) {
				optimalsS[t][0] = tOptTable[numIndex[0]][1];
				optimalsS[t][1] = tOptTable[numIndex[0] - 1][1] + tOptTable[numIndex[0] - 1][2];
				optimalsS[t][2] = tOptTable[numIndex[1]][1];
				optimalsS[t][3] = tOptTable[numIndex[1] - 1][1] + tOptTable[numIndex[1] - 1][2];
				optimalsS[t][4] = tOptTable[numIndex[2]][1];
				optimalsS[t][5] = tOptTable[numIndex[2] - 1][1] + tOptTable[numIndex[2] - 1][2];
			}
			else {
				int sIndex3 = numIndex[numIndex.length - 1];
				optimalsS[t][4] = tOptTable[sIndex3][1];
				optimalsS[t][5] = tOptTable[sIndex3 - 1][1] + tOptTable[sIndex3 - 1][2];
				int sIndex2 = numIndex[numIndex.length - 2];
				optimalsS[t][2] = tOptTable[sIndex2][1];
				optimalsS[t][3] = tOptTable[sIndex2 - 1][1] + tOptTable[sIndex2 - 1][2];
				int sIndex1 = numIndex[numIndex.length - 3];
				optimalsS[t][0] = tOptTable[sIndex1][1];
				optimalsS[t][1] = minSquare(optimalsS[t][0], sIndex1, tOptTable, t);
			}
		}
		System.out.println(Arrays.deepToString(optimalsS));
		return optimalsS;
	}
	
	
	double minSquare(double lb, int upIndex, double[][] tOptTable, int t) {
		int lowIndex = 0; double realS;
		for (int i = 0; i < tOptTable.length; i++)
			if (tOptTable[i][2] != maxOrderQuantity[t-1]) {
				lowIndex = i;
				break;
			}
		
		try{
			IloCplex cplex = new IloCplex();
			double ub = maxState;
			IloNumVar x  = cplex.numVar(lb, ub);
			realS = tOptTable[lowIndex][1] + tOptTable[lowIndex][2];
			IloNumExpr obj = cplex.prod(cplex.diff(x, realS), cplex.diff(x, realS));
			for (int i = lowIndex + 1; i <= upIndex; i++) {
				if (tOptTable[i][2] != maxOrderQuantity[t-1]) {
					realS = tOptTable[i][1] + tOptTable[i][2];
					obj = cplex.sum(obj, cplex.prod(cplex.diff(x, realS), cplex.diff(x, realS)));
				}
			}
			cplex.addMinimize(obj);
			if (cplex.solve()) {
				return cplex.getValue(x);
			}
		}
		catch (IloException e) {
	         System.err.println("Concert exception '" + e + "' caught");     
	    }
		return 0;
	}
	
	public static void writeToFile(String fileName, String str){
		File results = new File(fileName);
		try {
			FileOutputStream fos = new FileOutputStream(results, true);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			osw.write(str+"\n");
			osw.close();
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static String getHeadersString(){
		return "K,v,h,I0,pai,Qmax, DemandPatt, OpValue, Time(sec), simValue, error";
	}
	
	
	public static void main(String[] args) {
		writeToFile("./"+StructureDPforThreesS.class.getSimpleName() + "_results.csv", getHeadersString());
		
		double[][] demands = {{30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30},
				{6.6,9.3,11.1,12.9,16.8,21.6,24,26.4,31.5,33.9,36.3,40.8,44.4,47.1,48.3,50.1,50.1,48.9,47.1,44.4},
				{45.9,48.6,49.8,49.8,48.6,45.9,42.3,37.8,35.4,33,30.3,27.9,25.5,23.1,20.7,18.3,14.4,10.8,8.1,6.3},
				{36.3,30,23.7,21,23.7,30,36.3,39,36.3,30,23.7,21,23.7,30,36.3,30.9,24.3,21.3,26.4,33},
				{47.1,30,12.9,6,12.9,30,47.1,54,47.1,30,12.9,6,12.9,30,47.1,29.7,15,7.5,11.4,29.7},
				{62.7,27.3,9.9,23.7,1,22.8,32.7,34.5,67.2,6.6,14.4,40.8,3.9,63.3,25.5,45,53.1,24.6,10.2,49.8},
				{1,15.3,45.6,140.1,80.4,146.7,133.8,74.4,84.3,108.9,46.5,87.9,66,27.9,32.1,88.5,161.7,36,32.4,49.8},
				{14.1,24.3,70.8,118.2,49.2,86.1,152.4,117.3,226.2,208.2,78.3,58.5,96,33.3,57.3,115.5,17.7,134.7,127.8,179.7},
				{13.2,34.8,79.2,43.2,43.8,59.4,22.2,54.9,61.2,34.2,49.5,95.4,35.7,144.6,160.2,103.8,150.6,86.4,122.7,63.6},
				{14.7,56.4,19.2,83.7,135.9,67.2,66.9,155.1,87.3,164.1,193.8,67.2,64.5,132,34.8,131.4,132.9,35.7,172.5,152.1},
		};
		
		
		double[] K = {2000,1000,500};
		double[] v = {2,5,10};
		double[] pai = {15,10,5};
		double[] capacity = {3, 5, 1, 4, 3, 2, 1, 6, 4, 3, 2, 3, 4, 5, 6, 1, 2, 4, 3, 1};
		
		double truncationQuantile = 0.999;  
		double stepSize = 1; 
		double minState = -300;
		double maxState = 500;
		double holdingCost = 1;
		
		for (int iK = 0; iK < K.length; iK++) {
	    	  for (int iv = 0; iv < v.length; iv++) {
	    		  for (int ipai = 0; ipai < pai.length; ipai++) {
	    			  for (int idemand = 9; idemand < demands.length; idemand++) 
	    				  for ( int icapacity = 0; icapacity < capacity.length; icapacity++){	      

	   double[] meanDemand = demands[idemand];
	   double fixedOrderingCost = K[iK] ; 
	   double proportionalOrderingCost = v[iv]; 
	   double penaltyCost = pai[ipai];
	   int[] maxOrderQuantity = new int[meanDemand.length];
	   for (int i = 0; i < meanDemand.length; i++)
		   maxOrderQuantity[i] = (int) (Math.round(Arrays.stream(meanDemand).sum()/meanDemand.length)*capacity[icapacity]);
	   
	   StructureDPforThreesS inventory = new StructureDPforThreesS(fixedOrderingCost, proportionalOrderingCost, penaltyCost, 
				holdingCost, stepSize, minState, maxState, truncationQuantile,
				meanDemand, maxOrderQuantity);

		double initialInventory = 0;	      
		int period = 1;

		State initialState = inventory.new State(period, initialInventory);
		long currTime2=System.currentTimeMillis();
		double finalValue = inventory.f(initialState);
		System.out.println("final optimal expected value is: " + finalValue);
		
		double optQ = inventory.cacheActions.get(inventory.new State(period, initialInventory));
		System.out.println("optimal order quantity in the first priod is : " + optQ);	
		
		double time = (System.currentTimeMillis()-currTime2)/1000;
		System.out.println("running time is " + time + " s"); 
		
		Map<State, Double> cacheActions = inventory.cacheActions;
		Iterator<Map.Entry<State, Double>> iterator = cacheActions.entrySet().iterator();
		double[][] arr = new double[cacheActions.size()][3];
		int i = 0;
		while (iterator.hasNext()) {
			Map.Entry<State, Double> entry = iterator.next();
			//System.out.println(entry.getKey() + ",  Q = " + entry.getValue());
			arr[i++] =new double[]{entry.getKey().period, entry.getKey().initialInventory, entry.getValue()};
		} 		
		
		int sampleNum = 100000;
		double simFinalValue = inventory.simulateSamples(inventory.getsS(arr), inventory.generateSamples(meanDemand, sampleNum));
		System.out.println("final simulated expected value is: " + simFinalValue);
		System.out.printf("Optimality gap is: %.2f%%\n", (simFinalValue - finalValue)/finalValue*100);
		String out = fixedOrderingCost+",\t"+
                proportionalOrderingCost+",\t"+
                holdingCost+",\t"+
                initialInventory+",\t"+
                penaltyCost+",\t"+
                maxOrderQuantity+",\t"+
                (idemand + 1) +",\t"+
                finalValue +",\t"+
                time + ",\t" + 
                simFinalValue + ",\t" +
                (simFinalValue-finalValue)/finalValue;
         
        writeToFile("./"+StructureDPforThreesS.class.getSimpleName() + "_results.csv", out);
	}

	    		  }
	    	  }
		}
	}
}
