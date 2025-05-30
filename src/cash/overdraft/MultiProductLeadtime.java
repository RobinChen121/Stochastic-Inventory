package cash.overdraft;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.IntStream;

import sdp.cash.multiItem.Actions;
import sdp.cash.multiItem.CashRecursionMultiLead;
import sdp.cash.multiItem.CashStateMultiLead;
import sdp.cash.multiItem.Demands;
import sdp.cash.multiItem.GetPmfMulti;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.cash.multiItem.CashRecursionMulti;
import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
*@author: zhenchen
*@date: Apr 2, 2024, 7:24:49 PM
*@desp: multi product overdraft with leadtime 1 period;
*
* running very slow for normal demands; 2 periods gamma demands very fast;
* 3 periods or more for normal, gamma or Poison will be burdensome for DP, 3 hour no solution;
* 
* 3 periods for 3 value discrete distribution, final optimal cash  is 91.26875; final optimal cash  is 441.57499999999993 for overhead cost 0;
* final optimal cash  is 272.23749999999995 for overhead cost 50 (Q1 = 40, Q2 = 20);
* optimal order quantity in the first period is :  Q1 = 40, Q2 = 20
* running time is 137.0s;
 *
* double[][] values = {{20, 30, 40}, {10, 15, 20}};
* double[][] probs = {{0.25, 0.5, 0.25}, {0.25, 0.5, 0.25}};
 * final optimal cash  is 91.19499999999998
 * optimal order quantity in the first period is :  Q1 = 40, Q2 = 20
 * running time is 2863.0s
 *
 * when T = 2, final optimal cash  is -17.800000000000008
 * optimal order quantity in the first period is :  Q1 = 40, Q2 = 20
 * running time is 0.0s
 *
 * 3 periods:
 * double[][] values = {{10, 30}, {5, 15}};
 * double[][] probs = {{0.5, 0.5}, {0.5, 0.5}};
 * final optimal cash  is -76.56
 * optimal order quantity in the first period is :  Q1 = 30, Q2 = 15
 * running time is 1568.0s
* 
* double[] price = {5, 10};
		double[] variCost = {1, 2};  // lower margin vs higher margin
		double[] salValueUnit = Arrays.stream(variCost).map(a -> a*0.5).toArray();
		
		double iniCash = 0;  // initial cash
		int iniI1 = 0;  // initial inventory
		int iniI2 = 0;
		
		double r0 = 0; // deposite rate
		double r1 = 0.1;  // overdraft rate
		double r2 = 2; // penalty interest rate for overdraft exceeding the limit
		double limit = 500; // overdraft limit
		double interestFreeAmount = 0;
		
		double Qbound = 45; // maximum ordering quantity when having enough cash
		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 200;
		double minCashState = -500;
		double maxCashState = 5000;
		double discountFactor = 1;
		
		int T = 3; // horizon length
		double[] overheadCost = new double[T];
		Arrays.fill(overheadCost, 100); 
		double[] meanDemands = new double[] {30, 15};		
		double[][] demand = new double[2][T]; // higher average demand vs lower average demand
		double[] beta = {10, 1}; // lower variance vs higher variance
*
*/

public class MultiProductLeadtime {

	public static void main(String[] args) {
		double[] price = {5, 10};
		double[] variCost = {1, 2};  // lower margin vs higher margin
		double[] salValueUnit = Arrays.stream(variCost).map(a -> a*0.5).toArray();
		
		double iniCash = 0;  // initial cash
		int iniI1 = 0;  // initial inventory
		int iniI2 = 0;
		
		double r0 = 0; // deposit rate
		double r1 = 0.1;  // overdraft rate
		double r2 = 2; // penalty interest rate for overdraft exceeding the limit
		double limit = 500; // overdraft limit
		double interestFreeAmount = 0;
		
		double Qbound = 50; // maximum ordering quantity when having enough cash
		double truncationQuantile = 0.9999;
		int stepSize = 1;
		double minInventoryState = 0;
		double maxInventoryState = 200;
		double minCashState = -500;
		double maxCashState = 5000;
		double discountFactor = 1;
		
		// gamma distribution:mean demand is shape / beta and variance is shape / beta^2
		// beta = 1 / scale
		// shape = demand * beta
		// variance = demand / beta
		// gamma in ssj library: alpha is alpha, and lambda is beta(beta)
		int T = 3; // horizon length
		double[] overheadCost = new double[T];
		Arrays.fill(overheadCost, 100);
		double[] meanDemands = new double[] {10, 10};

		double[][] demand = new double[2][T]; // higher average demand vs lower average demand
		double[] beta = {10, 1}; // lower variance vs higher variance
		
		for (int t = 0; t < T; t++) {
			demand[0][t] = meanDemands[0];
			demand[1][t] = meanDemands[1];
		}
		
		int N = demand.length; // number of products
		// get demand possibilities for each period
//		Distribution[][] distributions =  new GammaDist[N][T];
//		Distribution[][] distributions =  new PoissonDist[N][T];
//		Distribution[][] distributions =  new NormalDist[N][T];
		double[][] values = {{10, 30}, {5, 15}};
		double[][] probs = {{0.5, 0.5}, {0.5, 0.5}};
		Distribution[][] distributions =  new DiscreteDistribution[N][T];
//		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
//		.mapToObj(i -> new DiscreteDistribution(values[i], probs[i], values[i].length)) // can be changed to other distributions
//		.toArray(DiscreteDistribution[]::new);
		
		for (int i = 0; i < N; i++)
			for (int t = 0; t < T; t++) {
//				 distributions[i][t] = new GammaDist(demand[i][t]* beta[i], beta[i]);
//				 distributions[i][t]= new NormalDist(demand[i][t], 0.5 * demand[i][t]);
//				 distributions[i][t] = new PoissonDist(demand[i][t]);
				distributions[i][t] = new DiscreteDistribution(values[i], probs[i], values[i].length);
			}	
		GetPmfMulti PmfMulti = new GetPmfMulti(distributions, truncationQuantile, stepSize);
		
		// build action list for two items
		Function<CashStateMultiLead, ArrayList<Actions>> buildActionList = s -> {
			ArrayList<Actions> actions = new ArrayList<>();
			for (int i = 0; i < Qbound; i++)
				for (int j = 0; j < Qbound; j++) {
					Actions thisAction = new Actions(i, j);
					actions.add(thisAction);				
				}
			return actions;
		};
		
		
		// Immediate Value Function	      
		ImmediateValueFunction<CashStateMultiLead, Actions, Demands, Double> immediateValue
		= (IniState, Actions, RandomDemands) -> {
			double action1 = Actions.getFirstAction();
			double action2 = Actions.getSecondAction();
			double demand1 = RandomDemands.getFirstDemand();
			double demand2 = RandomDemands.getSecondDemand();
			double preQ1 = IniState.getPreQ1();
			double preQ2 = IniState.getPreQ2();
			double endInventory1 = Math.max(0, IniState.getIniInventory1() + preQ1 - demand1);
			double endInventory2 = Math.max(0, IniState.getIniInventory2() + preQ2 - demand2);
//			double revenue1 = price[0] * (IniState.getIniInventory1() + preQ1 - endInventory1);
//			double revenue2 = price[1] * (IniState.getIniInventory2() + preQ2 - endInventory2);
			double revenue1 = price[0] * Math.min(demand1, IniState.getIniInventory1() + preQ1);
			double revenue2 = price[1] * Math.min(IniState.getIniInventory2() + preQ2 , demand2);
			double revenue = revenue1 + revenue2;
			double orderingCost1 = variCost[0] * action1;
			double orderingCost2 = variCost[1] * action2;
			double orderingCosts = orderingCost1 + orderingCost2;
			double salValue = 0;
			if (IniState.getPeriod() == T) {
				salValue = salValueUnit[0] * endInventory1 + salValueUnit[1] * endInventory2;
			}			
			int t = IniState.getPeriod() - 1;
			double cashBalanceBefore = IniState.getIniCash()- orderingCosts - overheadCost[t];// whether plus revenue in this time point
			double interest = 0;
			if (cashBalanceBefore >= 0)
				interest = -r0 * cashBalanceBefore;
			else if(cashBalanceBefore >= -interestFreeAmount)
				interest = 0;
			else if (cashBalanceBefore >= -limit) 
				interest = r1 * (-cashBalanceBefore - interestFreeAmount);
			else 
				interest = r2 * (-cashBalanceBefore - limit) + r1 * (limit - interestFreeAmount);

			double cashBalanceAfter = cashBalanceBefore - interest + revenue + salValue;
			double cashIncrement = cashBalanceAfter - IniState.getIniCash();
			return cashIncrement;
		};
		
		
		// State Transition Function
		StateTransitionFunction<CashStateMultiLead, Actions, Demands, CashStateMultiLead> stateTransition = (IniState, Actions, RandomDemands) -> {
			double preQ1 = IniState.getPreQ1();
			double preQ2 = IniState.getPreQ2();
			double action1 = Actions.getFirstAction();
			double action2 = Actions.getSecondAction();
			double nextPreQ1 = action1;
			double nextPreQ2 = action2;
			double endInventory1 = IniState.getIniInventory1() + preQ1 - RandomDemands.getFirstDemand();
			endInventory1 = Math.max(0, endInventory1);
			double endInventory2 = IniState.getIniInventory2() + preQ2 - RandomDemands.getSecondDemand();
			endInventory2 = Math.max(0, endInventory2);
			double nextCash = IniState.getIniCash() + immediateValue.apply(IniState, Actions, RandomDemands);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			endInventory1 = endInventory1 > maxInventoryState ? maxInventoryState : endInventory1;
			endInventory2 = endInventory2 < minInventoryState ? minInventoryState : endInventory2;
//			nextCash = (int) nextCash;  //  rounding states to save computing time
			endInventory1 = (int) endInventory1;
			endInventory2 = (int) endInventory2;
			return new CashStateMultiLead(IniState.getPeriod() + 1, endInventory1, endInventory2, nextPreQ1, nextPreQ2, nextCash);
		};
		
		
		/*******************************************************************
		 * Solve
		 */
		CashRecursionMultiLead recursion = new CashRecursionMultiLead(discountFactor, PmfMulti, buildActionList,
				                             stateTransition, immediateValue, T);
		int period = 1;
		CashStateMultiLead iniState = new CashStateMultiLead(period, iniI1, iniI2, 0, 0, iniCash);
		long currTime = System.currentTimeMillis();
		double finalValue = iniCash + recursion.getExpectedValue(iniState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first period is :  Q1 = " + recursion.getAction(iniState).getFirstAction()
				                      + ", Q2 = " + recursion.getAction(iniState).getSecondAction());
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");	

	}

}


