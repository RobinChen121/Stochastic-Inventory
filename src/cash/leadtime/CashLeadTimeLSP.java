package cash.leadtime;

import java.util.stream.DoubleStream;

import umontreal.ssj.probdist.PoissonDist;

/** 
* @author chen zhen 
* @version 创建时间：2018年4月8日 下午2:14:56 
* @value 类说明: 包含提前期的随机动态规划
*/
public class CashLeadTimeLSP {
	double price;
	double fixedOrderingCost, proportionalOrderingCost; 
	int stepSize;
	int leadTime;
	double minCashState, maxCashState, minInventoryState, maxInventoryState;
	double truncationQuantile;
	double[] demands;
	double[][][] pmf;
	double overheadCost;
	
	public CashLeadTimeLSP(double price, double fixedOrderingCost, double proportionalOrderingCost, 
			int stepSize, int leadtime, double minInventoryState, double maxInventoryState, double minCashState, double maxCashState,
			 double truncationQuantile, double[] demands, double overheadCost) {
		this.price = price;
		this.fixedOrderingCost = fixedOrderingCost;
		this.proportionalOrderingCost = proportionalOrderingCost;
		this.stepSize = stepSize;
		this.leadTime = leadtime;
		this.minInventoryState = minInventoryState;
		this.maxInventoryState = maxInventoryState;
		this.minCashState = minCashState;
		this.maxCashState = maxCashState;
		this.truncationQuantile = truncationQuantile;
		this.demands = demands;
		this.pmf = getPmf(demands);
		this.overheadCost = overheadCost;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
		
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
    	
		double cashIncrement = revenue - fixedCost - variableCost - overheadCost;
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

	public static void main(String[] args) {
		

	}

}
