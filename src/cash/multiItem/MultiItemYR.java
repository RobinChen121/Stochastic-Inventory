/**
 * @date: Jul 6, 2020
 */
package cash.multiItem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;

import sdp.cash.multiItem.CashRecursionV;
import sdp.cash.multiItem.CashSimulationMultiXR;
import sdp.cash.multiItem.CashSimulationY;
import sdp.cash.multiItem.CashStateMulti;
import sdp.cash.multiItem.CashStateMultiYR;
import sdp.cash.multiItem.CashStateR;
import sdp.cash.multiItem.GetPmfMulti;
import sdp.inventory.FinalCash.BoundaryFuncton;
import sdp.inventory.ImmediateValue.ImmediateValueFunctionV;
import sdp.inventory.StateTransition.StateTransitionFunctionV;
import sdp.write.WriteToExcel;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 6, 2020
 * @Desc: recode the dynamic programming of two-product problem by new states: (y1, y2, R) and two functional
 *        equations: V(y1, y2, R) and Pi(y1, y2, R)
 *
 */
public class MultiItemYR {
	
	
	
	
	
	public static void main(String[] args) {
		double[] price = {4, 10};
		double[] variCost = {2, 4};  // higher margin vs lower margin
		double depositeRate = 0;
		double[] salPrice = {1, 1};
		
		double iniCash = 20;  // initial cash
		int iniInventory1 = 0;  // initial inventory
		int iniInventory2 = 0;
			
		// gamma distribution:mean demand is shape * scale and variance is shape * scale^2
		// shape = demand * scale
		// variance = demand / scale
		double[][] demand = {{5, 5, 5}, {8, 8, 8}}; // higher average demand vs lower average demand
		double[] rate = {2, 1}; // higher variance vs lower variance
			
		int T = demand[0].length; // horizon length
		int m = demand.length; // number of products
		
		double truncationQuantile = 0.9999; // may affect poisson results
		int stepSize = 1;
		double minCashState = -50;
		double maxCashState = 10000;
		int minInventoryState = 0;	
		int maxInventoryState = 200;
		int Qbound = 20;
		double discountFactor = 1;
		
		// get demand possibilities for each period
		Distribution[][] distributions =  new GammaDist[m][T];
		//Distribution[][] distributions =  new PoissonDist[m][T];
		//Distribution[][] distributions =  new NormalDist[m][T];
		for (int i = 0; i < m; i++)
			for (int t = 0; t < T; t++) {
				distributions[i][t] = new GammaDist(demand[i][t]* rate[i], rate[i]);
				//distributions[i][t] = new PoissonDist(demand[i][t]);
				//distributions[i][t]= new NormalDist(demand[i][t], 0.1 * demand[i][t]);
			}
		
		// build action list (y1, y2) for V(x1, x2, R)
		Function<CashStateMultiYR, ArrayList<double[]>> buildActionListPai = s -> {
			ArrayList<double[]> actions = new ArrayList<>();
			double Ybound = Qbound;
			for (double i = 0; i < Ybound; i++)
				for (double j = 0; j < Ybound; j++) {
					double[] thisActions = {i, j};
					actions.add(thisActions);	
				}
			return actions;
		};
		
		// build action list (y1, y2) for Pai(x1, x2, R)
		Function<CashStateMulti, ArrayList<double[]>> buildActionListV = s -> {
			ArrayList<double[]> actions = new ArrayList<>();
			int miny1 = (int) s.getIniInventory1();
			int miny2 = (int) s.getIniInventory2();
			double iniR = s.getIniCash() + variCost[0] * s.getIniInventory1() + variCost[1] * s.getIniInventory2();
			for (int i = miny1; i < miny1 + Qbound; i++)
				for (int j = miny2; j < miny2 + Qbound; j++) {				
					if (variCost[0] * i + variCost[1] * j < iniR + 0.1) {
						double[] thisActions = {i, j};
						actions.add(thisActions);
					}					
				}
			return actions;
		};

		// Immediate Value Function	for V(y1, y2, R), the increment for R 
		
		ImmediateValueFunctionV<CashStateMultiYR, double[], Double> immediateValue
		= (IniState, RandomDemands) -> {
			double demand1 = RandomDemands[0];
			double demand2 = RandomDemands[1];
			double endInventory1 = Math.max(0, IniState.getIniInventory1()  - demand1);
			double endInventory2 = Math.max(0, IniState.getIniInventory2()  - demand2);
			double revenue1 = price[0] * Math.min(IniState.getIniInventory1(), demand1);
			double revenue2 = price[1] * Math.min(IniState.getIniInventory2(), demand2);
			double revenue = revenue1 + revenue2;
			//double initialCash = IniState.getIniR() - variCost[0] * IniState.getIniInventory1() - variCost[1] * IniState.getIniInventory2();
			double orderingCostY1 = variCost[0] * IniState.getIniInventory1();
			double orderingCostY2 = variCost[1] * IniState.getIniInventory2();
			double orderingCostsY = orderingCostY1 + orderingCostY2;
			double salValue = 0;
			if (IniState.getPeriod() == T) {
				salValue = salPrice[0] * endInventory1 + salPrice[1] * endInventory2;
			}
			return revenue + (1 - depositeRate) * (IniState.getIniR() - orderingCostsY) - IniState.getIniR();
		};
	
		BoundaryFuncton<CashStateMulti, Double> boundFinalCash
		= (IniState) -> {
			return IniState.getIniCash() + salPrice[0] * IniState.getIniInventory1() + salPrice[1] * IniState.getIniInventory2();
		};
	
	// State Transition Function
	StateTransitionFunctionV<CashStateMultiYR, double[], CashStateMulti> stateTransition
	= (IniState, RandomDemands) -> {
		double endInventory1 = IniState.getIniInventory1() - RandomDemands[0];
		endInventory1 = Math.max(0, endInventory1);
		double endInventory2 = IniState.getIniInventory2()- RandomDemands[1];
		endInventory2 = Math.max(0, endInventory2);
		double revenue1 = price[0] * Math.min(IniState.getIniInventory1(), RandomDemands[0]);
		double revenue2 = price[1] * Math.min(IniState.getIniInventory2(), RandomDemands[1]);
		double nextW = revenue1 + revenue2 + (1 + depositeRate) * (IniState.getIniR() - variCost[0] * IniState.getIniInventory1()
									- variCost[1] * IniState.getIniInventory2());  // revise
		
		endInventory1 = Math.round(endInventory1);
		endInventory2 = Math.round(endInventory2);
		nextW = Math.round(nextW);
		nextW = nextW > maxCashState ? maxCashState : nextW;
		nextW = nextW < minCashState ? minCashState : nextW;
		endInventory1 = endInventory1 > maxInventoryState ? maxInventoryState : endInventory1;
		endInventory2 = endInventory2 < minInventoryState ? minInventoryState : endInventory2;
		
		return new CashStateMulti(IniState.getPeriod() + 1, endInventory1, endInventory2, nextW);
	};

	
	
	GetPmfMulti PmfMulti = new GetPmfMulti(distributions, truncationQuantile, stepSize);
	
	/*******************************************************************
	 * Solve
	 */
	CashRecursionV recursion = new CashRecursionV(discountFactor, PmfMulti, buildActionListV, buildActionListPai,
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
	
	CashStateR iniState2 = new CashStateR(period, iniCash);
	double[] optY = recursion.getYStar(iniState2);
	System.out.println("optimal order quantity y* in the first priod is : " + Arrays.toString(optY));
	double[][] optTable = recursion.getOptTableDetail();
	
	
	/*******************************************************************
	 * Simulate
	 * 
	 * 
	 */
	int sampleNum = 10000;	
	currTime = System.currentTimeMillis();
	CashSimulationY simulation = new CashSimulationY(sampleNum, distributions, discountFactor, 
			 recursion, stateTransition);
	double simFinalValue = simulation.simulateSDPGivenSamplNum(iniState, variCost);
	double gap = (simFinalValue - finalValue) / finalValue;
	System.out.printf("optimality gap for this policy is %.2f%%\n", gap * 100);
	time = (System.currentTimeMillis() - currTime) / 1000.0;
	System.out.println("running time is " + time + "s");
	
	
	WriteToExcel wr = new WriteToExcel();
	String fileName = "Pai_yStar" + ".xls";
	String headString =  "period" + "\t" + "x1" + "\t" + "x2" + "\t" + "w" + "\t" + 
	          "c1" + "\t" + "c2" + "\t" + "R" + "\t" + "y1*"+ "\t" + "y2*" + "\t" + 
			   "cashConstrained" + "\t" + "alpha" + "\t" + "y1"  + "\t" + "y2";
	wr.writeArrayToExcel(optTable, fileName, headString);
	
	
	
	
	

}
	
}
