package cash.multiItem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;


import sdp.cash.multiItem.CashRecursionV2;
import sdp.cash.multiItem.CashStateMulti;
import sdp.cash.multiItem.GetPmfMulti;
import sdp.inventory.FinalCash.BoundaryFuncton;
import sdp.inventory.StateTransition.StateTransitionFunction;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;

/**
 * @author chen
 * @email: 15011074486@163.com
 * @Date: 2021 Feb 25, 19:48:25
 * @Description: multi item cash constrained problem for states (x1, x2, w)
 * 
 */
public class MultiItemCashXW {

	public static void main(String[] args) {
		double[] price = {2, 10};
		double[] variCost = {1, 2};  // higher margin vs lower margin
		double depositebeta = 0;		
		
		double iniCash = 10;  // initial cash
		int iniInventory1 = 0;  // initial inventory
		int iniInventory2 = 0;
			
		// gamma distribution:mean demand is shape / beta and variance is shape / beta^2
		// beta = 1 / scale
		// shape = demand * beta
		// variance = demand / beta
		// gamma in ssj: alpha is alpha, and lambda is beta(beta)
		int T = 4; // horizon length
		double[] meanDemands = new double[] {10, 3};
		
		double[][] demand = new double[2][T]; // higher average demand vs lower average demand
		double[] beta = {10, 1}; // higher variance vs lower variance
		
		double d1 = meanDemands[0];
		double d2 = meanDemands[1];
		double v1 = variCost[0]; double v2 = variCost[1];
		double p1 = price[0]; double p2 = price[1];
		for (int t = 0; t < T; t++) {
			demand[0][t] = d1;
			demand[1][t] = d2;
		}			
		
		double[] salPrice = Arrays.stream(variCost).map(a -> a*0.5).toArray();
		int m = demand.length; // number of products		
		
		double truncationQuantile = 0.9999; // may affect poisson results
		int stepSize = 1;
		double minCashState = 0;
		double maxCashState = 10000;
		int minInventoryState = 0;	
		int maxInventoryState = 200;
		int Qbound = 40;
		double discountFactor = 1;
		
		// get demand possibilities for each period
		Distribution[][] distributions =  new GammaDist[m][T];
		//Distribution[][] distributions =  new PoissonDist[m][T];
		//Distribution[][] distributions =  new NormalDist[m][T];
		for (int i = 0; i < m; i++)
			for (int t = 0; t < T; t++) {
				distributions[i][t] = new GammaDist(demand[i][t]* beta[i], beta[i]);
				//distributions[i][t] = new PoissonDist(demand[i][t]);
				//distributions[i][t]= new NormalDist(demand[i][t], 0.1 * demand[i][t]);
			}
		GetPmfMulti PmfMulti = new GetPmfMulti(distributions, truncationQuantile, stepSize);
		
		// build action list (y1, y2) for pai(x1, x2, w)
		// CashStateMulti are states (x1, x2, w)
		Function<CashStateMulti, ArrayList<double[]>> buildActionListPai = s -> {
			ArrayList<double[]> actions = new ArrayList<>();
			double Ybound = Qbound;
			for (double i = 0; i < Ybound; i=i+1)
				for (double j = 0; j < Ybound; j=j+1) {
					double[] thisActions = {i, j};
					actions.add(thisActions);	
				}
			return actions;
		};
		
		// build action list (y1, y2) for V(x1, x2, w)
		Function<CashStateMulti, ArrayList<double[]>> buildActionListV = s -> {
			ArrayList<double[]> actions = new ArrayList<>();
			int miny1 = (int) s.getIniInventory1();
			int miny2 = (int) s.getIniInventory2();
			double iniR = s.getIniCash() + v1 * s.getIniInventory1() + v2 * s.getIniInventory2();
			for (double i = miny1; i < miny1 + Qbound; i = i + 1)
				for (double j = miny2; j < miny2 + Qbound; j = j + 1) {				
					if (v1 * i + v2 * j < iniR + 0.1) {
						double[] thisActions = {i, j};
						actions.add(thisActions);
					}					
				}
			return actions;
		};		

	
		BoundaryFuncton<CashStateMulti, Double> boundFinalCash
		= (IniState) -> {
			return IniState.getIniCash() + salPrice[0] * IniState.getIniInventory1() + salPrice[1] * IniState.getIniInventory2();
		};
		
		// State Transition Function
		StateTransitionFunction<CashStateMulti, double[], double[], CashStateMulti> stateTransition
		= (IniState, actions, RandomDemands) -> {
			double endInventory1 = actions[0] - RandomDemands[0];
			endInventory1 = Math.max(0, endInventory1);
			double endInventory2 = actions[1] - RandomDemands[1];
			endInventory2 = Math.max(0, endInventory2);
			double revenue1 = p1 * Math.min(actions[0], RandomDemands[0]);
			double revenue2 = p2 * Math.min(actions[1], RandomDemands[1]);
			double nextW = revenue1 + revenue2 + (1 + depositebeta) * (IniState.getIniCash() - v1 * actions[0]
										- v2 * actions[1] + v1 * IniState.getIniInventory1() + v2 *  IniState.getIniInventory2());  // revise
			
			endInventory1 = Math.round(endInventory1 * 10) / 10.0;  // rounding the states to one decimal
			endInventory2 = Math.round(endInventory2 * 10) / 10.0;
			nextW = Math.round(nextW * 10) / 10.0;
			
			nextW = nextW > maxCashState ? maxCashState : nextW;
			nextW = nextW < minCashState ? minCashState : nextW;
			endInventory1 = endInventory1 > maxInventoryState ? maxInventoryState : endInventory1;
			endInventory2 = endInventory2 < minInventoryState ? minInventoryState : endInventory2;
			
			return new CashStateMulti(IniState.getPeriod() + 1, endInventory1, endInventory2, nextW);
		};
		
		
		/*******************************************************************
		 * Solve
		 */
		CashRecursionV2 recursion = new CashRecursionV2(discountFactor, PmfMulti, buildActionListV, buildActionListPai,
				stateTransition, boundFinalCash, T, variCost);
		int period = 1;
		CashStateMulti iniState = new CashStateMulti(period, iniInventory1, iniInventory2, iniCash);
		long currTime = System.currentTimeMillis();
		double finalValue = recursion.getExpectedValueV(iniState);
		System.out.println("final optimal cash  is " + finalValue);
		System.out.println("optimal order quantity in the first priod is :  y1 = " + recursion.getAction(iniState)[0]
				                      + ", y2 = " + recursion.getAction(iniState)[1]);
		double time = (System.currentTimeMillis() - currTime) / 1000.0;
		System.out.println("running time is " + time + "s");
		
		double[] optY = recursion.getYStar(iniState);
		System.out.println("optimal order quantity y* in the first priod is : " + Arrays.toString(optY));		

		
		/*******************************************************************
		 * Simulate
		 * 
		 * basically, this simulation is testing for Theorem 1: 
		 * optimal ordering decisions depend on y*(R)
		 */
		int sampleNum = 10000;	
		currTime = System.currentTimeMillis();

	}

}
