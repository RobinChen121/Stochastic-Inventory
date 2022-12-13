package cash.multiItem;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;

import sdp.cash.RecursionG;
import sdp.cash.multiItem.CashRecursionV2;
import sdp.cash.multiItem.CashSimulationY;
import sdp.cash.multiItem.CashStateMulti;
import sdp.cash.multiItem.GetPmfMulti;
import sdp.inventory.FinalCash.BoundaryFuncton;
import sdp.inventory.GetPmf;
import sdp.inventory.StateTransition.StateTransitionFunction;
import sdp.write.ReadExcel;
import sdp.write.WriteToExcelTxt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.UniformIntDist;

/**
 * @author chen
 * @email: 15011074486@163.com
 * @Date: 2021 Mar 1, 10:36:23
 * @Description: TODO 
 * 
 */
public class MultiItemCashXWTesting {
	public static void main(String[] args) {
		double[] price = {2, 10};
		double[] variCost = {1, 2};  // higher margin vs lower margin
		double depositeRate = 0;		
		
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
		
		// read parameter settings from excel files
		ReadExcel re = new ReadExcel();
		double[][] paraSettings = re.readExcelXLSX("Numerical experiments-2021-02-06.xlsx", 2);	
		
		for (int runTime = 1; runTime < 2; runTime=runTime+1) {
			price = new double[] {paraSettings[runTime][2], paraSettings[runTime][8]};
			variCost = new double[] {paraSettings[runTime][1], paraSettings[runTime][7]};
			double[] beta = new double[] {paraSettings[runTime][6], paraSettings[runTime][12]};
			// a = new int[] {(int)paraSettings[runTime][5], (int)paraSettings[runTime][11]}; // integer for uniform
			meanDemands = new double[] {paraSettings[runTime][3], paraSettings[runTime][9]};
		
		double[][] demand = new double[2][T]; // higher average demand vs lower average demand
		
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
		double stepSize = 1;
		double minCashState = 0;
		double maxCashState = 10000;
		int minInventoryState = 0;	
		int maxInventoryState = 200;
		int Qbound = 20;
		double discountFactor = 1;
		
		// get demand possibilities for each period
		//Distribution[][] distributions =  new GammaDist[m][T];
		//Distribution[][] distributions =  new PoissonDist[m][T];
		//Distribution[][] distributions =  new NormalDist[m][T];
		Distribution[][] distributions =  new UniformIntDist[m][T];
		for (int i = 0; i < m; i++)
			for (int t = 0; t < T; t++) {
				//distributions[i][t] = new GammaDist(demand[i][t]* beta[i], beta[i]);
				distributions[i][t] = new UniformIntDist((int)(demand[i][t] * 0.2), (int)(demand[i][t] * 1.8));
				//distributions[i][t] = new PoissonDist(demand[i][t]);
				//distributions[i][t]= new NormalDist(demand[i][t], 0.1 * demand[i][t]);
				// distributions[i][t] = new UniformIntDist(a[i], b[i]);
			}
		GetPmfMulti PmfMulti = new GetPmfMulti(distributions, truncationQuantile, stepSize);
		
		// build action list (y1, y2) for pai(x1, x2, w)
		// CashStateMulti are states (x1, x2, w)
		Function<CashStateMulti, ArrayList<double[]>> buildActionListPai = s -> {
			ArrayList<double[]> actions = new ArrayList<>();
			double Ybound = Qbound;
			for (double i = 0; i < Ybound; i = i + 1)
				for (double j = 0; j < Ybound; j = j + 1) {
					double[] thisActions = {i, j};
					actions.add(thisActions);	
				}
			return actions;
		};
		
		// build action list (y1, y2) for V(x1, x2, w), no use in the recursion
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
		StateTransitionFunction<CashStateMulti, double[], double[], CashStateMulti> stateTransition  // revise
		= (IniState, actions, RandomDemands) -> {
			double x1 = IniState.getIniInventory1();
			double x2 = IniState.getIniInventory2();
			double y1 = actions[0];
			double y2 = actions[1];
			double endInventory1 =  y1 - RandomDemands[0];
			endInventory1 = Math.max(0, endInventory1);
			double endInventory2 = y2  - RandomDemands[1];
			endInventory2 = Math.max(0, endInventory2);

			double revenue1 = p1 * Math.min(y1, RandomDemands[0]);
			double revenue2 = p2 * Math.min(y2, RandomDemands[1]);
			double nextW = revenue1 + revenue2 + (1 + depositeRate) * (IniState.getIniCash() - v1 * (y1-x1)
										- v2 * (y2-x2));  // revise
			
			endInventory1 = Math.round(endInventory1 * 1) / 1;  // rounding the states to one decimal 10.0
			endInventory2 = Math.round(endInventory2 * 1) / 1;  // very slow when decimal
			nextW = Math.round(nextW * 1) / 1;
			
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
		CashSimulationY simulation = new CashSimulationY(sampleNum, distributions, discountFactor, 
				 recursion, stateTransition);
		double simFinalValue = simulation.simulateSDPGivenSamplNum2(iniState, variCost);
		double gap = (finalValue - simFinalValue) / finalValue;
		System.out.printf("optimality gap for this policy y* is %.2f%%\n", gap * 100);
		time = (System.currentTimeMillis() - currTime) / 1000.0;
		System.out.println("running time is " + time + "s");
		
		
		/*******************************************************************
		 * Compute a1* and a2*
		 * 
		 * and simulate their results to test Theorem 2
		 * 
		*/
		stepSize = 1; // can be changed
		double[][][] pmf1 = new GetPmf(distributions[0], truncationQuantile, stepSize).getpmf();
		Distribution[] distributions1 = distributions[0];
		double[][][] pmf2 = new GetPmf(distributions[1], truncationQuantile, stepSize).getpmf();
		Distribution[] distributions2 = distributions[1];
		RecursionG recursionG1 = new RecursionG(pmf1, distributions1, price[0], variCost[0], 0, salPrice[0]);
		RecursionG recursionG2 = new RecursionG(pmf2, distributions2, price[1], variCost[1], 0, salPrice[1]);
		double[] opta1 = recursionG1.getOptY();
		double[] opta2 = recursionG2.getOptY();
		System.out.println("---------------");
		System.out.println("a1* in each period:");
		DecimalFormat df = new DecimalFormat("0.00");
		Arrays.stream(opta1).forEach(e -> System.out.print(df.format(e) + " " ));
		System.out.println("");
		System.out.println("a2* in each period:");
		Arrays.stream(opta2).forEach(e -> System.out.print(df.format(e) + " " ));
//		double simFinalValue2 = simulation.simulateSDPGivenSamplNuma1a22(iniState, variCost, opta1, opta2);
//		double gap2 = (simFinalValue2 - finalValue) / finalValue;
//		System.out.printf("optimality gap for this policy a* is %.2f%%\n", gap2 * 100);
		System.out.println();
		System.out.println("*****************************************");
		
		double[] mean = new double[] {demand[0][0], demand[1][0]};
		double[] variance = new double[] {demand[0][0] / beta[0], demand[1][0] / beta[1]};
		double[][] optTable = recursion.getOptTableDetail(mean, variance, price, opta1, opta2);
		
		double[] gaps = new double[] {gap};
		WriteToExcelTxt wr = new WriteToExcelTxt();
		String fileName = "run" + (int) paraSettings[runTime][0] + ".xls";
		String headString =  
					"meanD1" + "\t" + "meanD2" + "\t" + "variance1" + "\t" + "variance2" + "\t" +
			         "period" + "\t" + "x1" + "\t" + "x2" + "\t" + "w" + "\t" + 
					"p1" + "\t" + "p2" + "\t" +
			          "c1" + "\t" + "c2" + "\t" + "R" + "\t" + "y1*"+ "\t" + "y2*" + "\t" + 
					   "cashSituation" + "\t" + "alpha" + "\t" + "yHead1"  + "\t" + "yHead2"  + "\t" + "a1*"  + "\t" + "a2*" +
					   "\t" + "Theorem1Gap";
			wr.writeArrayToExcel2(optTable, fileName, headString, gaps);
		
	}
		
	}
}
