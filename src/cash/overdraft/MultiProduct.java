package cash.overdraft;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;

import sdp.cash.multiItem.Actions;
import sdp.cash.multiItem.CashStateMulti;
import sdp.cash.multiItem.Demands;
import sdp.cash.multiItem.GetPmfMulti;
import sdp.inventory.ImmediateValue.ImmediateValueFunction;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;

/**
*@author: zhenchen
*@date: Oct 17, 2023, 8:42:51 AM
*@desp: TODO
*
*/

public class MultiProduct {
	public static void main(String[] args) {
		double[] price = {5, 10};
		double[] variCost = {1, 2};  // higher margin vs lower margin
		double depositeRate = 0;
		
		
		double iniCash = 0;  // initial cash
		int iniI1 = 0;  // initial inventory
		int iniI2 = 0;
		
		double r0 = 0.01;
		double r1 = 0;
		double r2 = 0.1;
		double r3 = 1; // penalty interest rate for overdraft exceeding the limit
		double limit = 100; // overdraft limit
		double interestFreeAmount = 20;
		
		// gamma distribution:mean demand is shape / beta and variance is shape / beta^2
		// beta = 1 / scale
		// shape = demand * beta
		// variance = demand / beta
		// gamma in ssj library: alpha is alpha, and lambda is beta(beta)
		int T = 4; // horizon length
		double[] overheadCost = new double[T];
		Arrays.fill(overheadCost, 50); 
		double[] meanDemands = new double[] {10, 5};		
		double[][] demand = new double[2][T]; // higher average demand vs lower average demand
		double[] beta = {10, 1}; // lower variance vs higher variance
		
		for (int t = 0; t < T; t++) {
			demand[0][t] = meanDemands[0];
			demand[1][t] = meanDemands[1];
		}
				
		double[] salValueUnit = Arrays.stream(variCost).map(a -> a*0.5).toArray();
		int N = demand.length; // number of products
			
		double truncationQuantile = 0.9999; // may affect poisson results
		int stepSize = 1;
		double minCashState = 0;
		double maxCashState = 10000;
		int minInventoryState = 0;	
		int maxInventoryState = 200;
		int Qbound = 50;
		double discountFactor = 1;
		
		// get demand possibilities for each period
		Distribution[][] distributions =  new GammaDist[N][T];
		//Distribution[][] distributions =  new PoissonDist[m][T];
		//Distribution[][] distributions =  new NormalDist[m][T];
		for (int i = 0; i < N; i++)
			for (int t = 0; t < T; t++) {
				distributions[i][t] = new GammaDist(demand[i][t]* beta[i], beta[i]);
				//distributions[i][t] = new PoissonDist(demand[i][t]);
				//distributions[i][t]= new NormalDist(demand[i][t], 0.1 * demand[i][t]);
			}
		
		GetPmfMulti PmfMulti = new GetPmfMulti(distributions, truncationQuantile, stepSize);
		
		// build action list for two items
		Function<CashStateMulti, ArrayList<Actions>> buildActionList = s -> {
			ArrayList<Actions> actions = new ArrayList<>();
			for (int i = 0; i < Qbound; i++)
				for (int j = 0; j < Qbound; j++) {
					Actions thisAction = new Actions(i, j);
					actions.add(thisAction);				
				}
			return actions;
		};
		
		// Immediate Value Function	      
		ImmediateValueFunction<CashStateMulti, Actions, Demands, Double> immediateValue
		= (IniState, Actions, RandomDemands) -> {
			double action1 = Actions.getFirstAction();
			double action2 = Actions.getSecondAction();
			double demand1 = RandomDemands.getFirstDemand();
			double demand2 = RandomDemands.getSecondDemand();
			double endInventory1 = Math.max(0, IniState.getIniInventory1() + action1 - demand1);
			double endInventory2 = Math.max(0, IniState.getIniInventory2() + action2 - demand2);
			double revenue1 = price[0] * (IniState.getIniInventory1() + action1 - endInventory1);
			double revenue2 = price[1] * (IniState.getIniInventory2() + action2 - endInventory2);
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
				interest = r2 * (-cashBalanceBefore - interestFreeAmount);
			else 
				interest = r3 * (-cashBalanceBefore - limit) + r2 * (limit - interestFreeAmount);
			double cashBalanceAfter = cashBalanceBefore - interest + revenue + salValue;
			double cashIncrement = cashBalanceAfter - IniState.getIniCash();
			return cashIncrement;
		};
		
		// State Transition Function
		StateTransitionFunction<CashStateMulti, Actions, Demands, CashStateMulti> stateTransition = (IniState, Actions, RandomDemands) -> {
			double endInventory1 = IniState.getIniInventory1() + Actions.getFirstAction() - RandomDemands.getFirstDemand();
			endInventory1 = Math.max(0, endInventory1);
			double endInventory2 = IniState.getIniInventory2() + Actions.getSecondAction() - RandomDemands.getSecondDemand();
			endInventory2 = Math.max(0, endInventory2);
			double nextCash = IniState.getIniCash() + immediateValue.apply(IniState, Actions, RandomDemands);
			nextCash = nextCash > maxCashState ? maxCashState : nextCash;
			nextCash = nextCash < minCashState ? minCashState : nextCash;
			endInventory1 = endInventory1 > maxInventoryState ? maxInventoryState : endInventory1;
			endInventory2 = endInventory2 < minInventoryState ? minInventoryState : endInventory2;
			nextCash = (int) nextCash;  // rounding states to save computing time
			endInventory1 = (int) endInventory1;
			endInventory2 = (int) endInventory2;
			return new CashStateMulti(IniState.getPeriod() + 1, endInventory1, endInventory2, nextCash);
		};
		
		
		/*******************************************************************
		 * Solve
		 */
		CashRecursionMulti recursion = new CashRecursionMulti(discountFactor, PmfMulti, buildActionList,
				                             stateTransition, immediateValue, T);
		int period = 1;
		CashStateMulti iniState = new CashStateMulti(period, iniI1, iniI2, iniCash);
		long currTime = System.currentTimeMillis();
		double finalValue = iniCash + recursion.getExpectedValue(iniState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first priod is :  Q1 = " + recursion.getAction(iniState).getFirstAction()
				                      + ", Q2 = " + recursion.getAction(iniState).getSecondAction());
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");		
		
	}

}


