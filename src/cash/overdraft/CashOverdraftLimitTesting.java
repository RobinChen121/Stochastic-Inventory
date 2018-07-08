package cash.overdraft;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.DoubleStream;

import cash.constraint.CashConstraintLSPTesting;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.PoissonDist;

/** 
* @author chen zhen 
* @version 创建时间：2018年4月23日 下午9:38:22 
* @value 类说明: 用于大规模算例测试
*/
public class CashOverdraftLimitTesting{
	double price;
	double fixedOrderingCost, proportionalOrderingCost; 
	double holdingCost;
	double stepSize;
	double minCashState, maxCashState, minInventoryState, maxInventoryState;
	double truncationQuantile;
	double[] demands;
	double[][][] pmf;
	double interestRate;
	int maxQuantity;
	
	
	public CashOverdraftLimitTesting(double price, double fixedOrderingCost, double proportionalOrderingCost, 
			double holdingCost, double stepSize, double minInventoryState, double maxInventoryState, 
			double minCashState, double maxCashState, double truncationQuantile, 
			double[] demands, double interestRate, int maxQuantity) {
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
		this.interestRate = interestRate;
		this.maxQuantity = maxQuantity;
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
			int maxOrderQuantity = (int) Math.max((iniCash - minCashState - fixedOrderingCost)/ proportionalOrderingCost, 0);
			maxOrderQuantity = maxOrderQuantity > maxQuantity ? maxQuantity : maxOrderQuantity;
			return DoubleStream.iterate(0, i -> i + stepSize).limit(maxOrderQuantity + 1).toArray();
		}
		
		@Override
		public int hashCode(){
			String hash = "";
			hash = hash + period + 1 + iniInventory + 2 + iniCash;
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
			return "period = " + period +", "+"iniInventory = " + iniInventory + ", iniCash = " + iniCash; 
		}
	}
	
	double immediateValue(State state, double action, double randomDemand) {	
		double demand = randomDemand;
		double revenue = price*Math.min(state.iniInventory + action, demand);
    	double fixedCost = action > 0 ? fixedOrderingCost : 0;
    	double variableCost = proportionalOrderingCost*action;
    	double inventoryLevel = state.iniInventory + action -demand;
    	double holdCosts = holdingCost*Math.max(inventoryLevel, 0);
    	
    	double cashBalanceBefore = state.iniCash + revenue - fixedCost - variableCost - holdCosts;
    	double interest = interestRate*Math.max(-state.iniCash + fixedCost + variableCost, 0);
    	double cashBalanceAfter = cashBalanceBefore - interest;
    	double immediateValue = cashBalanceAfter - state.iniCash;
    	return immediateValue;
    }
	
	State stateTransition(State state, double action, double randomDemand) {
		double nextInventory = Math.max(0, state.iniInventory + action - randomDemand);
    	nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
    	nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
    	
    	double nextCash = state.iniCash + immediateValue(state, action, randomDemand);
    	nextCash = nextCash > maxCashState ? maxCashState : nextCash;
    	nextCash = nextCash < minCashState ? minCashState : nextCash;
    	
    	nextCash = Math.round(nextCash*10)/10.0;
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
			double val = -Double.MAX_VALUE; // min_value 是非零非负的最小值
			double bestOrderQty = 0;
			for (int i = 0; i < feasibleActions.length; i++) {
				double orderQty = feasibleActions[i];
				double thisQValue = 0;
				for (int j = 0; j < dAndP.length; j++) {
					thisQValue += dAndP[j][1]*immediateValue(s, orderQty, dAndP[j][0]);
					if (s.period < pmf.length)
						thisQValue += dAndP[j][1]*f(stateTransition(s, orderQty, dAndP[j][0]));
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
	
	double simulateSamples(State iniState, Map<State, Double> cacheActions, double[][] samples) {
		double[] cashs = new double[samples.length];
		for (int i = 0; i < samples.length; i++) {
			double[] I = new double[samples[0].length];
			double[] cashBalance = new double[samples[0].length];
			double[] demand = new double[samples[0].length];
			double Q = 0, fixCost, revenue, beforeBalance, interest;
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
						//System.out.println("error;");
						break;
						}
					I[t] = Math.max(0, I[t-1] + Q - demand[t]);
					revenue = price*Math.min(Q + I[t-1], demand[t]); 
				}
				fixCost = Q > 0 ? fixedOrderingCost : 0;
				if ( t== 0) {
					beforeBalance = iniState.iniCash + revenue - fixCost - proportionalOrderingCost*Q - holdingCost*I[t];
					interest = interestRate*Math.max(-iniState.iniCash + fixCost + proportionalOrderingCost*Q, 0);
				}
				else {	
					beforeBalance = cashBalance[t-1] + revenue - fixCost - proportionalOrderingCost*Q - holdingCost*I[t];
					interest = interestRate*Math.max(-cashBalance[t-1] + fixCost + proportionalOrderingCost*Q, 0);
				}				 
				cashBalance[t] = beforeBalance - interest;
				cashBalance[t] = cashBalance[t] > maxCashState ? maxCashState : cashBalance[t];
				cashBalance[t] = cashBalance[t] < minCashState ? minCashState : cashBalance[t];
				I[t] = I[t] > maxInventoryState ? maxInventoryState : I[t];
				I[t] = I[t] < minInventoryState ? minInventoryState : I[t];
				int period = t + 1;
				cashBalance[t] = Math.round(cashBalance[t]*10)/10.0;
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
			double Q, fixCost, demand, revenue, interest;
			for (int t = 0; t < samples[0].length; t++)
			{
				demand = samples[i][t];
				if ( t== 0) {
					Q = optsBS[t][2];
					I[t] = Math.max(0, Q - demand);
					revenue = price*Math.min(Q, demand); 
					fixCost = Q > 0 ? fixedOrderingCost : 0;
					cashFlow[t] = iniCash + revenue - fixCost-  proportionalOrderingCost*Q - holdingCost*I[t];
					interest = interestRate*Math.max(-iniCash + fixCost + proportionalOrderingCost*Q, 0);
					cashFlow[t] -= interest;
				}
				else {
					double maxOrderQuantity = Math.max(0, (cashFlow[t-1] - minCashState- fixedOrderingCost)/proportionalOrderingCost);
					maxOrderQuantity = maxOrderQuantity > maxQuantity ? maxQuantity : maxOrderQuantity;
					if (I[t-1] < optsBS[t][0] && cashFlow[t-1] > optsBS[t][1]) 
						Q = Math.min(maxOrderQuantity, optsBS[t][2] - I[t-1]);
					else
						Q = 0;
					I[t] = Math.max(0, I[t-1] + Q - demand);
					revenue = price*Math.min(Q + I[t-1], demand); 
					fixCost = Q > 0 ? fixedOrderingCost : 0;
					cashFlow[t] = cashFlow[t-1] + revenue - fixCost- proportionalOrderingCost*Q - holdingCost*I[t];
					interest = interestRate*Math.max(-cashFlow[t-1] + fixCost + proportionalOrderingCost*Q, 0);
					cashFlow[t] -= interest;
					cashFlow[t] = Math.round(cashFlow[t]*10)/10.0;
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
			int recordsTime = 0, recordBTime = 0;
			double[][] tOptTable = Arrays.stream(optimalTable).filter(p -> p[0] == i).map(p -> Arrays.stream(p).toArray())
											.toArray(double[][] :: new);
			optimalsS[t][2] = 0;
			optimalsS[t][1] = -80; // 初始B
			double mark_s = 0;
			for (int j = tOptTable.length - 1; j >= 0; j--) {
				if (tOptTable[j][3] != 0) {
					if (mark_s == 0) {
						optimalsS[t][0] = tOptTable[j+1][2]; // 确定小 s
						mark_s=1;
					}
					if (tOptTable[j][2] + tOptTable[j][3] > optimalsS[t][2])
						optimalsS[t][2] = tOptTable[j][2] + tOptTable[j][3]; // 确定大S,找一个最大的
					recordsTime = 1;
				}
				if (tOptTable[j][3] == 0 && recordsTime == 1) {
					if (tOptTable[j][1] > optimalsS[t][1])
						optimalsS[t][1] = tOptTable[j][1];
				}		
			}
		}
		System.out.println(Arrays.deepToString(optimalsS));
		return optimalsS;
	}
	
	public static String getHeadersString(){
	      return "K,v,S,b,demand,I0,B0,minB,Time(s),EC,ECSIM,sBSSIM,SimGap,sBSGap";
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
	
	public static void main(String[] args) {
		writeToFile("./"+CashOverdraftLimitTesting.class.getSimpleName() + "_results.csv", getHeadersString());

		double[] K = {10,15};
		double[] v = {1,2};
		double[] S = {5,10}; // selling price

		double[] minCash = {-40,-80};
		double[] b = {0.1,0.2};

		int initialCapital[] = {0,20};

		double[][] demands = {{7,7,7,7,7,7},
				{2,3,4,5,6,7},
				{8,7,6,5,4,3},
				{5,6,7,8,7,6},//7,6,5,4,5,6
				{8,5,2,1,2,5},
				{8,4,1,3,1,3},
				{1,3,8,4,8,7},
				{1,4,7,3,5,8},
				{3,8,4,4,6,2},
				{3,1,5,8,4,4}
		};

		for (int iB = 0; iB < initialCapital.length; iB++) 
			for (int iv = 0; iv < v.length; iv++) 
				for (int iK = 0; iK < K.length; iK++) 
					for (int is = 0; is < S.length; is++)
						for (int idemand = 0; idemand < demands.length; idemand++) 
							for (int ib = 0; ib < b.length; ib++)
								for (int iC = 0; iC < minCash.length; iC++){ 
		
		 
		double fixedOrderingCost = K[iK];
		double proportionalOrderingCost = v[iv];
		double holdingCost = 0;
		double price = S[is];		
		
		double iniCash = initialCapital[iB]; 
		double interestRate = b[ib];
				
		double truncationQuantile = 0.999;  
		int stepSize = 1;
		
		double minInventoryState = 0;
		double maxInventoryState = 100;
		double minCashState = minCash[iC]; // 资金不可能减少
		double maxCashState = 600;	
		int maxQuantity = 60;
		
		double[] meanDemands = demands[idemand];
		int T = meanDemands.length;
		PoissonDist[] distributions = new PoissonDist[T];
		for (int i = 0; i < T; i++) {
			distributions[i] = new PoissonDist(meanDemands[i]);
		}				
		
		CashOverdraftLimitTesting inventory = new CashOverdraftLimitTesting(price, fixedOrderingCost, proportionalOrderingCost, 
				holdingCost, stepSize, minInventoryState, maxInventoryState, 
				minCashState, maxCashState, truncationQuantile, meanDemands, interestRate, maxQuantity);
		
		int iniInventory = 0;	   
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
			//System.out.println(entry.getKey() + ",  Q = " + entry.getValue());
			optTable[i++] =new double[]{entry.getKey().period, entry.getKey().iniCash, entry.getKey().iniInventory, entry.getValue()};
		}
		int sampleNum = 10000;	
		Sampling sampling =new Sampling();
		double[][] samples = sampling.generateLHSamples(distributions, sampleNum);
		double simFinalValue = inventory.simulateSamples(initialState, cacheActions, samples);
		System.out.println("final simulated expected value is: " + simFinalValue);
		
		
		double[][] sBS = inventory.getsBS(iniCash, optTable);
		double simsSFinalValue = inventory.simulatesBS(iniCash, sBS, samples);
		System.out.println("final simulated sS expected value is: " + simsSFinalValue);
		System.out.printf("Optimality gap is: %.2f%%\n", (finalValue - simsSFinalValue)/finalValue*100);
		System.out.printf("Optimality gap is: %.2f%%\n", (simFinalValue - simsSFinalValue)/simFinalValue*100);
		
		String out = fixedOrderingCost+",\t"+
                proportionalOrderingCost+",\t"+
				price+",\t"+
                interestRate +",\t"+
				(idemand + 1) +",\t"+
                iniInventory+",\t"+
                iniCash+",\t"+
                minCashState+",\t"+ 
                time + ",\t" + 
                finalValue +",\t"+
                simFinalValue + ",\t" +
                simsSFinalValue + ",\t" +
                (finalValue - simsSFinalValue)/finalValue + ",\t" +
                (simFinalValue - simsSFinalValue)/simFinalValue;
         
        writeToFile("./"+ CashOverdraftLimitTesting.class.getSimpleName() + "_results.csv", out);
	}
}
}
