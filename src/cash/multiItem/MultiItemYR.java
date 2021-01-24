/**
 * @date: Jul 6, 2020
 */
package cash.multiItem;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;

import sdp.cash.RecursionG;
import sdp.cash.multiItem.CashRecursionV;
import sdp.cash.multiItem.CashSimulationMultiXR;
import sdp.cash.multiItem.CashSimulationY;
import sdp.cash.multiItem.CashStateMulti;
import sdp.cash.multiItem.CashStateMultiYR;
import sdp.cash.multiItem.CashStateR;
import sdp.cash.multiItem.GetPmfMulti;
import sdp.inventory.GetPmf;
import sdp.inventory.FinalCash.BoundaryFuncton;
import sdp.inventory.ImmediateValue.ImmediateValueFunctionV;
import sdp.inventory.StateTransition.StateTransitionFunctionV;
import sdp.write.ReadExcel;
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
		
//		for (int index = 5; index <= 10; index++) {
//			price[1] = index;
		
		
		double truncationQuantile = 0.9999; // may affect poisson results
		int stepSize = 1;
		double minCashState = 0;
		double maxCashState = 10000;
		int minInventoryState = 0;	
		int maxInventoryState = 200;
		int Qbound = 22;
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
		
		// build action list (y1, y2) for V(x1, x2, R)
		Function<CashStateMultiYR, ArrayList<double[]>> buildActionListPai = s -> {
			ArrayList<double[]> actions = new ArrayList<>();
			double Ybound = Qbound;
			for (double i = 0; i < Ybound; i=i+1)
				for (double j = 0; j < Ybound; j=j+1) {
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
	StateTransitionFunctionV<CashStateMultiYR, double[], CashStateMulti> stateTransition
	= (IniState, RandomDemands) -> {
		double endInventory1 = IniState.getIniInventory1() - RandomDemands[0];
		endInventory1 = Math.max(0, endInventory1);
		double endInventory2 = IniState.getIniInventory2()- RandomDemands[1];
		endInventory2 = Math.max(0, endInventory2);
		double revenue1 = p1 * Math.min(IniState.getIniInventory1(), RandomDemands[0]);
		double revenue2 = p2 * Math.min(IniState.getIniInventory2(), RandomDemands[1]);
		double nextW = revenue1 + revenue2 + (1 + depositebeta) * (IniState.getIniR() - v1 * IniState.getIniInventory1()
									- v2 * IniState.getIniInventory2());  // revise
		
		endInventory1 = Math.round(endInventory1 * 10) / 10;
		endInventory2 = Math.round(endInventory2 * 10) / 10;
		nextW = Math.round(nextW * 10) / 10;
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
	double[] mean = new double[] {demand[0][0], demand[1][0]};
	double[] variance = new double[] {demand[0][0] / beta[0], demand[1][0] / beta[1]};
	
	
	/*******************************************************************
	 * Simulate
	 * 
	 * basically, this simulation is testing for Theorem 1: 
	 * optimal ordering decisions depend on y*(R)
	 */
	int sampleNum = 10000;	
	currTime = System.currentTimeMillis();
	CashSimulationY simulation = new CashSimulationY(sampleNum, distributions, discountFactor, 
			 recursion, stateTransition);
	double simFinalValue = simulation.simulateSDPGivenSamplNum(iniState, variCost);
	double gap = (simFinalValue - finalValue) / finalValue;
	System.out.printf("optimality gap for this policy y* is %.2f%%\n", gap * 100);
	time = (System.currentTimeMillis() - currTime) / 1000.0;
	System.out.println("running time is " + time + "s");
	
	
	/*******************************************************************
	 * Compute a1* and a2*
	 * 
	 * and simulate their results to test Theorem 2
	 * 
	*/
	double[][][] pmf1 = new GetPmf(distributions[0], truncationQuantile, stepSize).getpmf();
	Distribution[] distributions1 = distributions[0];
	double[][][] pmf2 = new GetPmf(distributions[1], truncationQuantile, stepSize).getpmf();
	Distribution[] distributions2 = distributions[1];
	RecursionG recursionG1 = new RecursionG(pmf1, distributions1, price[0], variCost[0], 0, salPrice[0]);
	RecursionG recursionG2 = new RecursionG(pmf2, distributions2, price[1], variCost[1], 0, salPrice[1]);
	double[] opta1 = recursionG1.getOptY();
	double[] opta2 = recursionG2.getOptY();
	System.out.println("a1* in each period:");
	DecimalFormat df = new DecimalFormat("0.00");
	Arrays.stream(opta1).forEach(e -> System.out.print(df.format(e) + " " ));
	System.out.println("");
	System.out.println("a2* in each period:");
	Arrays.stream(opta2).forEach(e -> System.out.print(df.format(e) + " " ));
	double simFinalValue2 = simulation.simulateSDPGivenSamplNuma1a2(iniState, variCost, opta1, opta2);
	double gap2 = (simFinalValue2 - finalValue) / finalValue;
	System.out.printf("optimality gap for this policy a* is %.2f%%\n", gap2 * 100);
	double[][] optTable = recursion.getOptTableDetail2(mean, variance, price, opta1, opta2);
	
	double[] gaps = new double[] {gap, gap2};
	WriteToExcel wr = new WriteToExcel();
	String fileName = "run" + ".xls";
	String headString =  
			"meanD1" + "\t" + "meanD2" + "\t" + "variance1" + "\t" + "variance2" + "\t" +
	         "period" + "\t" + "x1" + "\t" + "x2" + "\t" + "w" + "\t" + 
			"p1" + "\t" + "p2" + "\t" +
	          "c1" + "\t" + "c2" + "\t" + "R" + "\t" + "y1*"+ "\t" + "y2*" + "\t" + 
			   "cashSituation" + "\t" + "alpha" + "\t" + "yHead1"  + "\t" + "yHead2"  + "\t" + "a1*"  + "\t" + "a2*" +
			   "\t" + "Theorem1Gap" + "Theorem2Gap";
	wr.writeArrayToExcel2(optTable, fileName, headString, gaps);
	
//	System.out.println("alpha in the first period: " + optTable[0][10]);
//	System.out.println("*******************************");
	
		}	
	
	
	

}
	
//}
