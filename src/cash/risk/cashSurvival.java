/**
 * @date: Nov 21, 2020
 */
package cash.risk;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.cash.CashRecursion;
import sdp.cash.CashSimulation;
import sdp.cash.CashState;
import sdp.cash.RiskRecursion;
import sdp.cash.RiskSimulation;
import sdp.cash.RiskState;
import sdp.cash.CashRecursion.OptDirection;
import sdp.inventory.GetPmf;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.sampling.Sampling;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Nov 21, 2020
 * @Desc: stochastic dynamic programming model to maximize the survival probability, 
 *        the paper is Archibald & Betts (2002), management science.
 *        more than 5 periods will run very slow. 6 periods about 300s of running time.
 *        rolling length 3, 10 samples each period, running time is 450s.
 *		  rolling length 2, 20 samples each period, running time is 69s.
 */
public class cashSurvival {

	/**
	 * @param args
	 * @date: Nov 21, 2020, 6:01:10 PM 
	 */
	public static void main(String[] args) {
		
		WriteToExcelTxt wr = new WriteToExcelTxt();
		String fileName = "SurvivalDiffCash.xls";
		String headString =  "iniCash" + "\t" + "optQ" + "\t" + "survivalProb" + "\t" + "serviceLeval";		
						
		double[] meanDemand = {14,23,33,46,50};
		int T = meanDemand.length;
		
		double iniI = 0;
		double iniCash = 80;
		double fixOrderCost = 0;
		double[] price = new double[T];
		double[] variCost = new double[T];
		Arrays.fill(price, 4);	
		Arrays.fill(variCost, 1);	
		double depositeRate = 0;
		double salvageValue = 0.5;
		double holdingCost = 0;	
		
		double[] overheadCosts = new double[T];	
		Arrays.fill(overheadCosts, 100); // overhead costs
		double overheadRate = 0; // rate from revenue to pay overhead wages
		double maxOrderQuantity = 1000; // maximum ordering quantity
		
		double truncationQuantile = 0.99;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 1000;
		double minCashState = -500; // can affect results, should be smaller than minus fixedOrderCost
		double maxCashState = 5000;
		double penaltyCost = 0; // can also be overdraft rate; large penalty cost cause big gaps for simulation results, since may generate zero demand;
		double discountFactor = 1;		
		
		int sampleNumSim = 300; // number of scenarios for rolling SAA
		int sampleNumPeriod = 15;
		int rollingLength = 1;
		double meanDemandSum = Arrays.stream(meanDemand).sum();
		double rollingDemandSum = Arrays.stream(meanDemand).limit(rollingLength).sum();
		double portion = rollingDemandSum / meanDemandSum;
		double serviceRateRequired = 0.8;
		double serviceRate = 0.95;
		double rollingServiceRate = Math.pow(serviceRate, portion);
		
		// get demand possibilities for each period
		double sigmaCoe = 0.2;
		
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
//				.mapToObj(i -> new NormalDist(meanDemand[i], sigmaCoe * meanDemand[i])).toArray(Distribution[]::new);// can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i])).toArray(Distribution[]::new);
		
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();
		
//		for (int k = 0; k < 40; k++) {
//			iniCash = k * 10;

		// feasible actions
		Function<RiskState, double[]> getFeasibleAction = s -> {
			int t = s.getPeriod() - 1;
			double maxQ = Math.min(s.iniCash/variCost[t], maxOrderQuantity);
			if (s.getBankruptBefore() == true)
				maxQ = 0;
			maxQ = Math.max(maxQ, 0);
			return DoubleStream.iterate(0, i -> i + stepSize).limit((int) maxQ + 1).toArray();
		};	
		
		int period = 1;	
		RiskState initialState = new RiskState(period, iniI, iniCash, false);
		
		// immediate value
		ImmediateValueFunction<RiskState, Double, Double, Double> immediateValue = (state, action, randomDemand) -> {
			int t = state.getPeriod() - 1;
			double revenue = price[t] * Math.min(state.getIniInventory() + action, randomDemand);
			double fixedCost = action > 0 ? fixOrderCost : 0;
			double variableCost = variCost[t] * action;
			double deposite = (state.getIniCash() - fixedCost - variableCost) * (1 + depositeRate);
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdCosts = holdingCost * Math.max(inventoryLevel, 0);
			double cashIncrement = revenue + deposite - holdCosts - overheadCosts[t]
					- state.getIniCash();
			double salValue = state.getPeriod() == T ? salvageValue * Math.max(inventoryLevel, 0) : 0;
			cashIncrement += salValue;
			return cashIncrement;
		};
		
		// state transition function
		StateTransitionFunction<RiskState, Double, Double, RiskState> stateTransition = (state, action,
				randomDemand) -> {
			double nextInventory = Math.max(0, state.getIniInventory() + action - randomDemand);
			double nextCash = state.getIniCash() + immediateValue.apply(state, action, randomDemand);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			nextInventory = nextInventory > maxInventoryState ? maxInventoryState : nextInventory;
			nextInventory = nextInventory < minInventoryState ? minInventoryState : nextInventory;
			// cash is integer or not
			nextCash = Math.round(nextCash * 1) / 1; // the right should be a decimal
			boolean bankruptBefore = state.getBankruptBefore();
			if (nextCash < 0)
				bankruptBefore = true;
			return new RiskState(state.getPeriod() + 1, nextInventory, nextCash, bankruptBefore);
		};
		
		long currTime = System.currentTimeMillis();
		NumberFormat nf = NumberFormat.getPercentInstance();
		nf.setMinimumFractionDigits(5);		
		DecimalFormat df = new DecimalFormat("###, ###");
		
		/*******************************************************************
		 * Solve SDP
		 */
		RiskRecursion recursion = new RiskRecursion(pmf, getFeasibleAction, stateTransition, immediateValue);
		recursion.setTreeMapCacheAction();
		double finalValue = recursion.getSurvProb(initialState);
		System.out.println("survival probability for this initial state is: " + finalValue);
		double optQ = recursion.getAction(initialState);
		System.out.println("optimal order quantity in the first priod is : " + optQ);
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
//		double[][] optTable = recursion.getOptTable();
		System.out.println();
		
		/*******************************************************************
		 * Simulate the result
		 */
		
		int sampleNum = 1000;
		RiskSimulation simulation = new RiskSimulation(distributions, sampleNum, recursion); // no need to add overheadCost in this class
		double[] result = simulation.simulateLostSale(initialState, immediateValue);
		DecimalFormat df2 = new DecimalFormat("###, ###");
		System.out.println("\nfinal simulated survival probability in " + df2.format(sampleNum) + " samples is: " + nf.format(result[0]));
		System.out.println("\nfinal simulated lost sale rate " + " is: " + nf.format(result[1]));
		System.out.println("************************************************************");
		
		/*******************************************************************
		 * solve the problem by SDP when there is individual chance constraint approximation
		 */			
		// feasible actions 2
		Function<RiskState, double[]> getFeasibleAction2 = s -> {
			int t = s.getPeriod() - 1;
			double thisPeriodServRate = 1- (1 - serviceRateRequired) / T;
			double cashQ = Math.max(0, s.iniCash/variCost[t]);
			double minQ = Math.ceil(distributions[t].inverseF(thisPeriodServRate)); 
			minQ = minQ > cashQ ? cashQ : minQ;
			double maxQ = Math.min(cashQ, maxOrderQuantity);
//			if (s.getBankruptBefore() == true)
//				maxQ = 0;
			minQ = Math.ceil(distributions[t].inverseF(thisPeriodServRate));
			maxQ = cashQ;
			if (maxQ < minQ)
				return new double[]{cashQ};
			else {
				int K = (int) maxQ - (int) minQ + 1;
				double[] actions = new double[K];
				for (int i = 0; i < K; i++)
					actions[i] = (int) minQ + i;
				return actions;
			}
			
		};
		
		/*******************************************************************
		 * Solve
		 */		
		recursion = new RiskRecursion(pmf, getFeasibleAction2, stateTransition, immediateValue);
		period = 1;		
		initialState = new RiskState(period, iniI, iniCash, false);
		currTime = System.currentTimeMillis();
		recursion.setTreeMapCacheAction();
		double SDPLbSurvival;
		SDPLbSurvival= recursion.getSurvProb(initialState);
		System.out.println("**********************************************");
		System.out.println("result of SDP with service rate constraint is: ");
		System.out.println("survival probability for this initial state is: " + nf.format(SDPLbSurvival));
		System.out.println("optimal order quantity in the first priod is : " + recursion.getAction(initialState));
		time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		/*******************************************************************
		 * Simulating sdp results
		 */		
		
		sampleNum = 1000;
		simulation = new RiskSimulation(distributions, sampleNum, recursion); // no need to add overheadCost in this class
		result = simulation.simulateLostSale(initialState, immediateValue);
		System.out.println("final simulated survival probability in " + df.format(sampleNum) + " samples is: " + nf.format(result[0]));
		System.out.println("final simulated lost sale rate " + " is: " + nf.format(result[1]));
		double serviceSDPLB = 1 - result[1];
		System.out.println("final simulated service rate " + " is: " + nf.format(serviceSDPLB));
		
		/*******************************************************************
		 * rolling horizon approach
		 * 
		 */
		int[] sampleNumSims = new int[T]; // sample number in each period, the number of samples in the first period can have big influence		
		int[] sampleNumSimsRolling = new int[T];
		Arrays.fill(sampleNumSimsRolling, sampleNumPeriod);
	    currTime = System.currentTimeMillis();
	    System.out.println("**********************************************");
	    // generate samples
	 	Sampling sampling = new Sampling();
	 	double[][] scenarios = sampling.generateLHSamples(distributions, sampleNumSims);
	    scenarios = sampling.generateLHSamples(distributions, sampleNumSimsRolling);
	    
	   
	    
	    double[] resultSim; 
	    RiskSimulation simulation1 = new RiskSimulation(distributions, sampleNumSim, immediateValue, stateTransition); // no need to add overheadCost in this class
	    
	    resultSim = simulation1.rollingHoirzon(rollingLength, initialState, rollingServiceRate, sampleNumSimsRolling, price, variCost, overheadCosts, salvageValue, holdingCost, scenarios, sampleNumSim);
	    double time1 = (System.currentTimeMillis() - currTime) / 1000.00;	    
	    System.out.println("after rolling horizon for length " + rollingLength +", " + "total horizon length is " + T + ", result is: ");
	    System.out.println("running time is " + time1 + "s");
	    double rollingObj = resultSim[0];
	    System.out.println("final simulated survival probability of rolling SAA in " + sampleNumSim + " samples is: " + nf.format(rollingObj));
	    double sigma2 = Math.sqrt(resultSim[1]*(1 - resultSim[1])/sampleNumSim);
		double error2  = 1.96*sigma2;
		double serviceRateRolling = 1 - resultSim[1];
		double optQ1Rolling = resultSim[2];
		System.out.println("the service rate for simulated SAA rolling horizon is " + nf.format(serviceRateRolling));
		System.out.println("the optimal ordering Q in the 1st period of the SAA rolling horizon is " + optQ1Rolling);
		
		
//		double[] out = new double[]{iniCash, optQ, finalValue, serviceRate};
//		wr.writeToExcelAppend(out, fileName);
//		 
//		}
	}

}







