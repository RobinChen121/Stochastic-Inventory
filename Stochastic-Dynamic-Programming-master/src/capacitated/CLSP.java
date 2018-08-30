package capacitated;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import umontreal.ssj.probdist.PoissonDist;

public class CLSP {
	double[][][] pmf;
	
	public CLSP(double[][][] pmf) {
		this.pmf = pmf;
	}

	class State{
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
	
	Function<State, double[]> actionGenerator;
	
    
	@FunctionalInterface
	interface StateTransitionFunction <S, A, R, S2>{
		public S2 apply (S s, A a, R r);
	}
	StateTransitionFunction <State, Double, Double, State>  stateTransition;
	
	@FunctionalInterface
	interface ImmediateValueFunction <S, A, R, V>{
		public V apply (S s, A a, R r);
	}
	ImmediateValueFunction<State, Double, Double, Double> immediateValue;
		
	
	Comparator<State> keyComparator = (o1, o2) -> o1.period > o2.period ? 1 : 
		o1.period == o2.period ? o1.initialInventory > o2.initialInventory ? 1 : 
			(o1.initialInventory == o2.initialInventory ? 0 : -1) : -1;
	SortedMap<State, Double> cacheActions = new TreeMap<>(keyComparator);
	SortedMap<State, Double> cacheValues = new TreeMap<>(keyComparator);
	
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
	
	int[] levelNum(double[][] optTable, int maxOrderQuantity) {
		ArrayList<Integer> indexArr = new ArrayList<>();
		boolean mark = false;
		for (int j = 0; j < optTable.length; j++) {
			if (optTable[j][2] < maxOrderQuantity && !mark) {
				mark = true;
			}
			else if (optTable[j][2] == maxOrderQuantity && mark) {
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
					
	
	public static void main(String[] args) {
		
		  
	      double initialInventory = 0; 
	      double[] meanDemand = {50, 60, 40};
	      
	      double truncationQuantile = 0.99;  // 置信度稍微改一下，泊松分布特殊
	      double stepSize = 1; 
	      double minState = -200;
	      double maxState = 300;
	      int T = meanDemand.length;

	      double fixedOrderingCost = 500; 
	      double proportionalOrderingCost = 0; 
	      double penaltyCost = 10;
	      double holdingCost = 2;
	      int maxOrderQuantity = 100;
	      
	      PoissonDist[] distributions = IntStream.iterate(0, i -> i + 1)
	                                              .limit(T)
	                                              .mapToObj(i -> new PoissonDist(meanDemand[i]))
	                                              .toArray(PoissonDist[]::new); // replace for loop
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

	      CLSP inventory = new CLSP(pmf);
	      
	      inventory.actionGenerator=s->{
	      	  return DoubleStream.iterate(0,i->i+stepSize).limit(maxOrderQuantity+1).toArray();
	        };
	      
	       inventory.stateTransition= (state, action, randomDemand) ->{
	    	   double nextInventory = state.initialInventory + action - randomDemand;
	    	   nextInventory = nextInventory > maxState ? maxState : nextInventory;
	    	   nextInventory = nextInventory < minState ? minState : nextInventory;
	    	   return inventory.new State(state.period+1, nextInventory);
	       };
	     

	       inventory.immediateValue= (state, action, randomDemand) ->
	      {
	    	  double fixedCost = action > 0 ? fixedOrderingCost : 0;
	    	  double variableCost = proportionalOrderingCost*action;
	    	  double inventoryLevel = state.initialInventory + action -randomDemand;
	    	  double holdingCosts = holdingCost*Math.max(inventoryLevel, 0);
	    	  double penaltyCosts = penaltyCost*Math.max(-inventoryLevel, 0);
	    	  double totalCosts = fixedCost + variableCost +holdingCosts + penaltyCosts;
	    	  return totalCosts;
	      };
	     
	      int period = 1;
	      State initialState = inventory.new State(period, initialInventory);	      
	      long currTime2=System.currentTimeMillis();
	      
	      double finalValue = inventory.f(initialState);
	      System.out.println("final optimal expected value is: " + finalValue);
	      
	      double optQ = inventory.cacheActions.get(inventory.new State(period, initialInventory));
	      System.out.println("optimal order quantity in the first priod is : " + optQ);
	      
	      Map<State, Double> cacheActions = inventory.cacheActions;
	      Iterator<Map.Entry<State, Double>> iterator = cacheActions.entrySet().iterator();
	      double[][] optTable = new double[cacheActions.size()][3];
	      int i = 0;
	      while (iterator.hasNext()) {
	    	  Map.Entry<State, Double> entry = iterator.next();
	    	  optTable[i++] =new double[]{entry.getKey().period, entry.getKey().initialInventory, entry.getValue()};
	      }
	      
	      
	      double time = (System.currentTimeMillis()-currTime2)/1000;
	      System.out.println("running time is " + time + " s"); 
	       
	}
}
