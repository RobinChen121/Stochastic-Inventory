package milp;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.IntUnaryOperator;
import java.util.stream.IntStream;


import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.CplexI;
import ilog.cplex.IloCplex;
import sdp.inventory.State;
import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

/** 
* @author: Zhen Chen
* @email: okchen321@163.com
* @date: 2018��10��19��, ����8:47:59 
* @email: okchen321@163.com
* @copyright: MIT licence
* @description: a heuristic method to obtain values of s, C, S for cash constrained 
*               stochastic lot sizing problem by solving the approximate deterministic 
*               problem (a mixed integer programming model);
*               tests show its performance is avg 1% gap.
*               
*               this class fits best for poisson distribution.
*               
*               
* @note: this class need cplex.jar              
* 
*/

public class MipCashConstraint {
	double iniInventory;
	double iniCash;		
	double fixOrderCost;
	double variCost;
	double holdingCost;
	double price;
	double salvageValue;
	double overheadCost;
	Distribution[] distributions;
	public Map<State, Double> cacheC1Values = new TreeMap<>(); // record C1 values for different initial inventory x

	
	public MipCashConstraint(double iniInventory, double iniCash, double fixOrderCost, double variCost, double holdingCost, 
			double price, double salvageValue, Distribution[] distributions, double overheadCost) {
		this.iniInventory = iniInventory;
		this.iniCash = iniCash;
		this.fixOrderCost = fixOrderCost;
		this.variCost = variCost;
		this.holdingCost = holdingCost;
		this.price = price;
		this.salvageValue = salvageValue;
		this.distributions = distributions;
		this.overheadCost = overheadCost;
		Comparator<State> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? 0 : -1 : -1;
		this.cacheC1Values = new TreeMap<>(keyComparator);
	}
	
	public MipCashConstraint(double iniInventory, double iniCash, double fixOrderCost, double variCost, double holdingCost, 
			double price, double salvageValue, Distribution[] distributions) {
		this.iniInventory = iniInventory;
		this.iniCash = iniCash;
		this.fixOrderCost = fixOrderCost;
		this.variCost = variCost;
		this.holdingCost = holdingCost;
		this.price = price;
		this.salvageValue = salvageValue;
		this.distributions = distributions;
		Comparator<State> keyComparator = (o1, o2) -> o1.getPeriod() > o2.getPeriod() ? 1 : 
			o1.getPeriod() == o2.getPeriod() ? o1.getIniInventory() > o2.getIniInventory() ? 1 : 
				o1.getIniInventory() == o2.getIniInventory() ? 0 : -1 : -1;
		this.cacheC1Values = new TreeMap<>(keyComparator);
	}
	
	/**
	 *  this mip model take lost sale as a decision variable
	 * @return s, C, S
	 */
	public double[][] findsCSLostSale() {		
		double[] varx = null;
		double[] vary;
		double[] varw;
		double[] varI = null;
		double[] varB = null;;
		
		try {
			int T = distributions.length;
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); // no cplex logging information
			
			// parameter values in array
			double[] S = new double[T];
			double[] h = new double[T];
			double[] v = new double[T];
			double[] p = new double[T];
			Arrays.fill(S, fixOrderCost);
			Arrays.fill(h, holdingCost);
			Arrays.fill(v, variCost);
			Arrays.fill(p, price);
					
			// decision variables
			IloIntVar[] x = cplex.boolVarArray(T);  // whether ordering in period t
			IloNumVar[] y = cplex.numVarArray(T, 0.0, Double.MAX_VALUE);  // how much to order in period t
			IloNumVar[] w = cplex.numVarArray(T, 0.0, Double.MAX_VALUE);  // how much lost sales in period t
			IloNumVar[] I = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); 
			IloNumVar[] B = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // end-of-period cash in each period
			
			// objective function
			IloNumExpr finalCash = cplex.sum(cplex.prod(salvageValue, I[T - 1]), B[T - 1]);
			cplex.addMaximize(finalCash);
			
			// constraints
			// inventory equality: I_t = I_{t-1} + y_t - (d_t - w_t)
			// cash flow: B_{t} = B_{t - 1} + p_t(d_t - w_t) - hI_t - vy_t - Sx_t
			// relationship between x_t and y_t
			// cash constraint
			
			IloNumExpr realDemand = cplex.numExpr();
			IloNumExpr tRevenue = cplex.numExpr();
			IloNumExpr tHoldCost= cplex.numExpr();
			IloNumExpr tPurchaseCost = cplex.numExpr();
			IloNumExpr tFixCost = cplex.numExpr();
			IloNumExpr tTotalCost = cplex.numExpr();
			IloNumExpr tTotalOrderCost = cplex.numExpr();
			
			for (int i = 0; i < T; i++) {
				realDemand = cplex.diff(distributions[i].getMean(), w[i]);
				tRevenue = cplex.prod(p[i], realDemand);
				tHoldCost = cplex.prod(h[i], I[i]);
				tPurchaseCost = cplex.prod(v[i], y[i]);
				tFixCost = cplex.prod(S[i], x[i]);
				tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
				tTotalCost = cplex.sum(tHoldCost, tTotalOrderCost);
				
				if (i == 0) {
					cplex.addEq(cplex.diff(I[i], iniInventory), cplex.diff(y[i], realDemand));
					cplex.addEq(cplex.diff(B[i], iniCash), cplex.diff(tRevenue, tTotalCost));
					cplex.addLe(tTotalOrderCost, iniCash);
				}
				else {
					cplex.addEq(cplex.diff(I[i], I[i - 1]), cplex.diff(cplex.sum(y[i], w[i]), distributions[i].getMean()));
					cplex.addEq(cplex.diff(B[i], B[i - 1]), cplex.diff(tRevenue, tTotalCost));
					cplex.addLe(tTotalOrderCost, B[i - 1]);
				}
				cplex.addLe(y[i], cplex.prod(x[i], 10000));
			}
			
			if (cplex.solve()) {				
				cplex.output().println("Solution status = " + cplex.getStatus());
				cplex.output().println("Solution value = " + cplex.getObjValue());
				varx = cplex.getValues(x);
				vary = cplex.getValues(y);
				varw = cplex.getValues(w);
				varI = cplex.getValues(I);
				varB = cplex.getValues(B);
				System.out.println("x = ");
				System.out.println(Arrays.toString(varx));
				System.out.println("y = ");
				System.out.println(Arrays.toString(vary));
				System.out.println("w = ");
				System.out.println(Arrays.toString(varw));
				System.out.println("I = ");
				System.out.println(Arrays.toString(varI));
				System.out.println("B = ");
				System.out.println(Arrays.toString(varB));
				
			}
			cplex.end();
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		
		return heuristicFindsCS(varx, varI, varB);
	}
	
	/**
	 * this mip model is for lost sale, lost sale quantity is not a decision variable;
	 * order-up-to level is a decision variable in this model.
	 * 
	 * @return s, C, S
	 */
	public double[][] findsCS() {		
		double[] varx = null;
		double[] vars;
		double[] varI = null;
		double[] varB = null;
		
		try {
			int T = distributions.length;
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); // no cplex logging information
			
			// parameter values in array
			double[] K = new double[T];
			double[] h = new double[T];
			double[] v = new double[T];
			double[] p = new double[T];
			Arrays.fill(K, fixOrderCost);
			Arrays.fill(h, holdingCost);
			Arrays.fill(v, variCost);
			Arrays.fill(p, price);
					
			// decision variables
			IloIntVar[] x = cplex.boolVarArray(T);  // whether ordering in period t
			// IloIntVar[] z = new IloIntVar[T]; 
			IloNumVar[] s = cplex.numVarArray(T, 0.0, Double.MAX_VALUE);  // order-up-to level in period t
			IloNumVar[] I = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); 
			IloNumVar[] B = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // end-of-period cash in each period
			
			// objective function
			IloNumExpr finalCash = cplex.sum(cplex.prod(salvageValue, I[T - 1]), B[T - 1]);
			cplex.addMaximize(finalCash);
			
			// constraints
			// cash constraint: B_t >= Kx_t + v(s_t - I_{t-1})
			// cash flow: B_{t} = B_{t - 1} + p_t(s_t - I_t) - hI_t - v(s_t - I_{t-1}) - Kx_t
			// s_t >= I_{t-1}
			// s_t <= I_{t} + d_{t}
			// s_t-I_{t-1} <= x_t*10000
			
			IloNumExpr realDemand = cplex.numExpr();
			IloNumExpr tRevenue = cplex.numExpr();
			IloNumExpr tHoldCost= cplex.numExpr();
			IloNumExpr tPurchaseCost = cplex.numExpr();
			IloNumExpr tFixCost = cplex.numExpr();
			IloNumExpr tTotalCost = cplex.numExpr();
			IloNumExpr tTotalOrderCost = cplex.numExpr();
			
			for (int i = 0; i < T; i++) {	
				realDemand = cplex.diff(s[i], I[i]);
				tRevenue = cplex.prod(p[i], realDemand);
				tHoldCost = cplex.prod(h[i], I[i]);				
				tFixCost = cplex.prod(K[i], x[i]);
							
				if (i == 0) {
					tPurchaseCost = cplex.prod(v[i], cplex.diff(s[i], iniInventory));
					tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
					tTotalCost = cplex.sum(tHoldCost, cplex.sum(tTotalOrderCost, overheadCost));
					cplex.addLe(iniInventory, s[i]);
					cplex.addEq(cplex.diff(B[i], iniCash), cplex.diff(tRevenue, tTotalCost));
					cplex.addLe(cplex.sum(overheadCost, tTotalOrderCost), iniCash);
					cplex.addGe(distributions[i].getMean(), cplex.diff(s[i], I[i]));
					cplex.addLe(cplex.diff(s[i], iniInventory), cplex.prod(x[i], 10000));
				}
				else {
					tPurchaseCost = cplex.prod(v[i], cplex.diff(s[i], I[i - 1]));
					tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
					tTotalCost = cplex.sum(tHoldCost, cplex.sum(tTotalOrderCost, overheadCost));
					cplex.addLe(I[i - 1], s[i]); // s_t >= I_{t-1}
					cplex.addEq(cplex.diff(B[i], B[i - 1]), cplex.diff(tRevenue, tTotalCost)); // B_{t} = B_{t - 1} + p_t(s_t - I_t) - hI_t - v(s_t - I_{t-1}) - Kx_t
					cplex.addLe(cplex.sum(overheadCost, tTotalOrderCost), B[i - 1]); // B_t >= Kx_t + v(s_t - I_{t-1})
					cplex.addGe(distributions[i].getMean(), cplex.diff(s[i], I[i])); // s_t-I_t <= d_t
					cplex.addLe(cplex.diff(s[i], I[i-1]), cplex.prod(x[i], 10000)); // s_t-I_{t-1} <= x_t*10000
				}				
			}
			
			if (cplex.solve()) {				
				cplex.output().println("Solution status = " + cplex.getStatus());
				cplex.output().println("Solution value = " + cplex.getObjValue());
				varx = cplex.getValues(x);
				vars = cplex.getValues(s);
				varI = cplex.getValues(I);
				varB = cplex.getValues(B);
				System.out.println("x = ");
				System.out.println(Arrays.toString(varx));
				System.out.println("order-up-to level = ");
				System.out.println(Arrays.toString(vars));
				System.out.println("I = ");
				System.out.println(Arrays.toString(varI));
				System.out.println("B = ");
				System.out.println(Arrays.toString(varB));
				
			}
			cplex.end();
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}

		return heuristicFindsCS(varx, varI, varB);
	}
	
	/**
	 * this mip model is for lost sale, lost sale quantity is not a decision variable;
	 * order-up-to level is a decision variable in this model.
	 * 
	 * revise the above model, but it seems the results of the two linear models are very similar
	 * 
	 * @return s, C, S
	 */
	public double[][] findsCSNew() {		
		double[] varx = null;
		double[] vars;
		double[] varI = null;
		double[] varB = null;
		
		try {
			int T = distributions.length;
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); // no cplex logging information
			
			// parameter values in array
			double[] K = new double[T];
			double[] h = new double[T];
			double[] v = new double[T];
			double[] p = new double[T];
			Arrays.fill(K, fixOrderCost);
			Arrays.fill(h, holdingCost);
			Arrays.fill(v, variCost);
			Arrays.fill(p, price);
					
			// decision variables
			IloIntVar[] x = cplex.boolVarArray(T);  // whether ordering in period t
			IloIntVar[] delta = cplex.boolVarArray(T);  // whether it is lost sale in period t
			IloNumVar[] s = cplex.numVarArray(T, 0.0, Double.MAX_VALUE);  // order-up-to level in period t
			IloNumVar[] I = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); 
			IloNumVar[] B = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // end-of-period cash in each period
			
			// objective function
			IloNumExpr finalCash = cplex.sum(cplex.prod(salvageValue, I[T - 1]), B[T - 1]);
			cplex.addMaximize(finalCash);
			
			// constraints
			// cash constraint: B_t >= Kx_t + v(s_t - I_{t-1})
			// cash flow: B_{t} = B_{t - 1} + p_t(s_t - I_t) - hI_t - v(s_t - I_{t-1}) - Kx_t
			// s_t >= I_{t-1}
			// s_t - I_{t-1} <= delta*M
			// may should also add: s_t-d_t <= delta_t*M
			// s_t - d_t <= I_{t} + (1-delta_t)M
			// s_t - d_t >= I_{t} - (1-delta_t)M
			// I_t <= delta_t M
			
			IloNumExpr realDemand = cplex.numExpr();
			IloNumExpr tRevenue = cplex.numExpr();
			IloNumExpr tHoldCost= cplex.numExpr();
			IloNumExpr tPurchaseCost = cplex.numExpr();
			IloNumExpr tFixCost = cplex.numExpr();
			IloNumExpr tTotalCost = cplex.numExpr();
			IloNumExpr tTotalOrderCost = cplex.numExpr();
			
			for (int i = 0; i < T; i++) {	
				realDemand = cplex.diff(s[i], I[i]);
				tRevenue = cplex.prod(p[i], realDemand);
				tHoldCost = cplex.prod(h[i], I[i]);				
				tFixCost = cplex.prod(K[i], x[i]);
							
				if (i == 0) {
					tPurchaseCost = cplex.prod(v[i], cplex.diff(s[i], iniInventory));
					tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
					tTotalCost = cplex.sum(tHoldCost, cplex.sum(tTotalOrderCost, overheadCost));
					cplex.addLe(iniInventory, s[i]);
					cplex.addEq(cplex.diff(B[i], iniCash), cplex.diff(tRevenue, tTotalCost));
					cplex.addLe(cplex.sum(overheadCost, tTotalOrderCost), iniCash);
					cplex.addGe(distributions[i].getMean(), cplex.diff(s[i], I[i]));
					cplex.addLe(cplex.diff(s[i], iniInventory), cplex.prod(x[i], 10000));
					cplex.addGe(cplex.diff(s[i], iniInventory), 0);
				}
				else {
					tPurchaseCost = cplex.prod(v[i], cplex.diff(s[i], I[i - 1]));
					tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
					tTotalCost = cplex.sum(tHoldCost, cplex.sum(tTotalOrderCost, overheadCost));
					cplex.addLe(I[i - 1], s[i]); // s_t >= I_{t-1}
					cplex.addEq(cplex.diff(B[i], B[i - 1]), cplex.diff(tRevenue, tTotalCost)); // B_{t} = B_{t - 1} + p_t(s_t - I_t) - hI_t - v(s_t - I_{t-1}) - Kx_t
					cplex.addLe(cplex.sum(overheadCost, tTotalOrderCost), B[i - 1]); // B_t >= Kx_t + v(s_t - I_{t-1})
					cplex.addGe(distributions[i].getMean(), cplex.diff(s[i], I[i])); // s_t-I_t <= d_t
					cplex.addLe(cplex.diff(s[i], I[i-1]), cplex.prod(x[i], 10000)); // s_t-I_{t-1} <= x_t*10000
					cplex.addGe(cplex.diff(s[i], I[i-1]), 0); // s_t - I_{t-1} <= delta*M
					cplex.addLe(I[i], cplex.prod(delta[i], 10000)); //I_t <= delta_t * M
					cplex.addLe(cplex.diff(s[i],distributions[i].getMean()), cplex.sum(I[i], cplex.prod(10000, cplex.diff(1, delta[i])))); // s_t - d_t <= I_{t} + (1-delta_t)M
					cplex.addGe(cplex.diff(s[i],distributions[i].getMean()), cplex.diff(I[i], cplex.prod(10000, cplex.diff(1, delta[i])))); // s_t - d_t >= I_{t} - (1-delta_t)M	
				}				
			}
			
			if (cplex.solve()) {				
				cplex.output().println("Solution status = " + cplex.getStatus());
				cplex.output().println("Solution value = " + cplex.getObjValue());
				varx = cplex.getValues(x);
				vars = cplex.getValues(s);
				varI = cplex.getValues(I);
				varB = cplex.getValues(B);
				System.out.println("x = ");
				System.out.println(Arrays.toString(varx));
				System.out.println("order-up-to level = ");
				System.out.println(Arrays.toString(vars));
				System.out.println("I = ");
				System.out.println(Arrays.toString(varI));
				System.out.println("B = ");
				System.out.println(Arrays.toString(varB));
				
			}
			cplex.end();
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}

		return heuristicFindsCS(varx, varI, varB);
	}
	
	/**
	 * compute L(y) in a single period
	 * @param y : order-up-to level y,
	 * @param t : period t + 1
	 * @param distribution: distribution at period t
	 */
	double Ly(double y, int t, Distribution distribution) {		
		double meanI = 0;
		for (int i = 0; i < y; i++)
			meanI += (y - i) * (distribution.cdf(i + 0.5) - distribution.cdf(i - 0.5));
		
		double Ly = 0;
		if (t == distributions.length - 1)
			Ly = (price - variCost) * y- (price + holdingCost - salvageValue) * meanI;
		else
			Ly = (price - variCost) * y- (price + holdingCost) * meanI;

		return Ly;
	}

	
	int nextIndex(double[] x, int t) {
		int i;

		for (i = t + 1; i < x.length; i++) {
			if (x[i] > 0.1) {
				//i =  i - 1; 
				i = x[0]== 0 && t== 0 ? i : i-1;
				break;
			}
//			double S = distributions[i].inverseF((price - variCost) / (holdingCost  + price));
//			if (Ly(S, t, distributions[i]) < fixOrderCost) {
//				i =  i - 1; 
//				break;
//			}
		}
		
		return Math.min(i, x.length - 1);
	}
	
	/**
	 *  this mip model consider penalty cost for negative cash
	 * @return s, C, S
	 */
	public double[][] findsCSNewPenalty(double penaltyCost) {		
		double[] varx = null;
		double[] vars;
		double[] varI = null;
		double[] varB = null;
		
		try {
			int T = distributions.length;
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); // no cplex logging information
			
			// parameter values in array
			double[] K = new double[T];
			double[] h = new double[T];
			double[] v = new double[T];
			double[] p = new double[T];
			Arrays.fill(K, fixOrderCost);
			Arrays.fill(h, holdingCost);
			Arrays.fill(v, variCost);
			Arrays.fill(p, price);
					
			// decision variables
			IloIntVar[] x = cplex.boolVarArray(T);  // whether ordering in period t
			IloIntVar[] delta = cplex.boolVarArray(T);  // whether it is lost sale in period t
			IloIntVar[] delta2 = cplex.boolVarArray(T);  // whether it is negative cash in the end of period t
			IloNumVar[] s = cplex.numVarArray(T, 0.0, Double.MAX_VALUE);  // order-up-to level in period t
			IloNumVar[] I = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); 
			IloNumVar[] B = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // end-of-period cash in each period
			
			// objective function
			IloNumExpr finalCash = cplex.sum(cplex.prod(salvageValue, I[T - 1]), B[T - 1]);
			cplex.addMaximize(finalCash);
			
			// constraints
			// cash constraint: B_t >= Kx_t + v(s_t - I_{t-1})
			// cash flow: B_{t} = B_{t - 1} + p_t(s_t - I_t) - hI_t - v(s_t - I_{t-1}) - Kx_t
			// s_t >= I_{t-1}
			// s_t - I_{t-1} <= delta*M
			// s_t - d_t <= I_{t} + (1-delta_t)M
			// s_t - d_t >= I_{t} - (1-delta_t)M
			// I_t <= delta_t M
			
			IloNumExpr realDemand = cplex.numExpr();
			IloNumExpr tRevenue = cplex.numExpr();
			IloNumExpr tHoldCost= cplex.numExpr();
			IloNumExpr tPurchaseCost = cplex.numExpr();
			IloNumExpr tFixCost = cplex.numExpr();
			IloNumExpr tTotalCost = cplex.numExpr();
			IloNumExpr tTotalOrderCost = cplex.numExpr();
			
			for (int i = 0; i < T; i++) {	
				realDemand = cplex.diff(s[i], I[i]);
				tRevenue = cplex.prod(p[i], realDemand);
				tHoldCost = cplex.prod(h[i], I[i]);				
				tFixCost = cplex.prod(K[i], x[i]);
				
				cplex.addGe(cplex.diff(B[i], B[i]), cplex.prod(penaltyCost, cplex.diff(B[i], cplex.prod(1000000, cplex.diff(1, delta[i])))));
				cplex.addLe(cplex.diff(B[i], B[i]), cplex.prod(penaltyCost, cplex.diff(B[i], cplex.prod(1000000, cplex.diff(1, delta[i])))));
				cplex.addLe(B[i], cplex.prod(10000000, cplex.diff(1, delta[i])));
				
				if (i == 0) {
					tPurchaseCost = cplex.prod(v[i], cplex.diff(s[i], iniInventory));
					tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
					tTotalCost = cplex.sum(tHoldCost, cplex.sum(tTotalOrderCost, overheadCost));
					cplex.addLe(iniInventory, s[i]);
					cplex.addEq(cplex.diff(B[i], iniCash), cplex.diff(tRevenue, tTotalCost));
					cplex.addLe(cplex.sum(overheadCost, tTotalOrderCost), iniCash);
					cplex.addGe(distributions[i].getMean(), cplex.diff(s[i], I[i]));
					cplex.addLe(cplex.diff(s[i], iniInventory), cplex.prod(x[i], 10000));
					cplex.addGe(cplex.diff(s[i], iniInventory), 0);
				}
				else {
					tPurchaseCost = cplex.prod(v[i], cplex.diff(s[i], I[i - 1]));
					tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
					tTotalCost = cplex.sum(tHoldCost, cplex.sum(tTotalOrderCost, overheadCost));
					cplex.addLe(I[i - 1], s[i]); // s_t >= I_{t-1}
					cplex.addEq(cplex.diff(B[i], B[i - 1]), cplex.diff(tRevenue, tTotalCost)); // B_{t} = B_{t - 1} + p_t(s_t - I_t) - hI_t - v(s_t - I_{t-1}) - Kx_t
					cplex.addLe(cplex.sum(overheadCost, tTotalOrderCost), B[i - 1]); // B_t >= Kx_t + v(s_t - I_{t-1})
					cplex.addGe(distributions[i].getMean(), cplex.diff(s[i], I[i])); // s_t-I_t <= d_t
					cplex.addLe(cplex.diff(s[i], I[i-1]), cplex.prod(x[i], 10000)); // s_t-I_{t-1} <= x_t*10000
					cplex.addGe(cplex.diff(s[i], I[i-1]), 0); // s_t - I_{t-1} <= delta*M
					cplex.addLe(I[i], cplex.prod(delta[i], 10000)); //I_t <= delta_t * M
					cplex.addLe(cplex.diff(s[i],distributions[i].getMean()), cplex.sum(I[i], cplex.prod(10000, cplex.diff(1, delta[i])))); // s_t - d_t <= I_{t} + (1-delta_t)M
					cplex.addGe(cplex.diff(s[i],distributions[i].getMean()), cplex.diff(I[i], cplex.prod(10000, cplex.diff(1, delta[i])))); // s_t - d_t >= I_{t} - (1-delta_t)M	
				}
				
			}
			
			if (cplex.solve()) {				
				cplex.output().println("Solution status = " + cplex.getStatus());
				cplex.output().println("Solution value = " + cplex.getObjValue());
				varx = cplex.getValues(x);
				vars = cplex.getValues(s);
				varI = cplex.getValues(I);
				varB = cplex.getValues(B);
				System.out.println("x = ");
				System.out.println(Arrays.toString(varx));
				System.out.println("order-up-to level = ");
				System.out.println(Arrays.toString(vars));
				System.out.println("I = ");
				System.out.println(Arrays.toString(varI));
				System.out.println("B = ");
				System.out.println(Arrays.toString(varB));
				
			}
			cplex.end();
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}

		return heuristicFindsCS(varx, varI, varB);
	}
	
	
	/**
	 * this mip model is for lost sale, lost sale quantity is not a decision variable;
	 * order-up-to level is a decision variable in this model.
	 * 
	 * I^+ is represented by piecewise model
	 * 
	 * revise the above model, but it seems the results of the two linear models are very similar
	 * 
	 * @return s, C, S
	 */
	public double[][] findsCSPieceWise() {
		int T = distributions.length;
		double[] varz = null;
		double[] varS;
		double[] varx = null;
		double[] varxplus = null;
		double[] varR = null;
		int M = 1000000;
		
		try {
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); // no cplex logging information
			
			// parameter values in array
			double[] K = new double[T];
			double[] h = new double[T];
			double[] v = new double[T];
			double[] p = new double[T];
			Arrays.fill(K, fixOrderCost);
			Arrays.fill(h, holdingCost);
			Arrays.fill(v, variCost);
			Arrays.fill(p, price);
					
			// decision variables
			IloIntVar[] z = cplex.boolVarArray(T);  // whether ordering in period t
			IloNumVar[][] P = new IloNumVar[T][T];  // whether the ordering cycle is from period i to period j
			for (int i = 0; i < P.length; i++)
				for (int j = 0; j < P.length; j++)
					P[i][j] = cplex.boolVar();	
			IloIntVar[] delta = cplex.boolVarArray(T);  // whether it is lost sale in period t
			IloNumVar[] S = cplex.numVarArray(T, 0.0, Double.MAX_VALUE);  // order-up-to level in period t
			IloNumVar[] x = cplex.numVarArray(T, -Double.MAX_VALUE, Double.MAX_VALUE);  // inventory
			IloNumVar[] xplus = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); 
			IloNumVar[] xminus = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // may be useless
			IloNumVar[] R = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // end-of-period cash in each period
			
			
			// objective function
			IloNumExpr finalCash = cplex.sum(cplex.prod(salvageValue, xplus[T - 1]), R[T - 1]);
			cplex.addMaximize(finalCash);
			
			IloNumExpr realDemand = cplex.numExpr();
			IloNumExpr tRevenue = cplex.numExpr();
			IloNumExpr tHoldCost= cplex.numExpr();
			IloNumExpr tPurchaseCost = cplex.numExpr();
			IloNumExpr tFixCost = cplex.numExpr();
			IloNumExpr tTotalCost = cplex.numExpr();
			IloNumExpr tTotalOrderCost = cplex.numExpr();
			
			int segmentNum = 5;
			double[] probs = new double[segmentNum];
			double segmentNumFloat = (float) segmentNum;
			double prob = (float) 1.0 /segmentNumFloat;
			Arrays.fill(probs, prob);
			int nbSamples = 10000;
			
			double[][] sumLambdas = new double[T][T]; // sum lambdas from period i to period j
			for (int i  =  0;  i < T; i++)
				for (int j = 0; j < T; j++) {
					if (i > j)
						sumLambdas[i][j] = 0;
					else {
						for (int k = i; k <= j; k++)
							sumLambdas[i][j] += distributions[k].getMean();
					}
				}
			
			
			for (int t = 0; t < T; t++) {	
				realDemand = cplex.diff(S[t], xplus[t]);
				tRevenue = cplex.prod(p[t], realDemand);
				tHoldCost = cplex.prod(h[t], xplus[t]);				
				tFixCost = cplex.prod(K[t], z[t]);
				

				// sum Pjt == 1
				IloLinearNumExpr sumPjt;
				IloLinearNumExpr sumzjt;
				IloLinearNumExpr sumzjt2;
				sumPjt = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) 
					sumPjt.addTerm(1, P[j][t]);
				cplex.addEq(sumPjt, 1); // upper triangle
				for (int j = t + 1; j < T; j++)
					cplex.addEq(P[j][t], 0);  // other Pjt = 0, or else cannot output values
				
				// Pjt >= z_j - sum_{k = j+1}^{t}z_k
				// sum_{j=1}{t}z_j == 0 => P[0][t] == 1, this constraints are important for the extra piecewise constraints
				sumzjt2 = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) {
					sumzjt = cplex.linearNumExpr();
					for (int k = j + 1; k <= t; k++)
						sumzjt.addTerm(z[k], 1);
					cplex.addGe(P[j][t], cplex.diff(z[j], sumzjt));
					sumzjt2.addTerm(z[j], 1);					
				}
				cplex.addGe(sumzjt2, cplex.diff(1, P[0][t]));
				
				// piecewise constraints
				// in each period t, get its conditional expectations for S_{jt}, j = 1, 2, ..., t
				long[] seed = {1,2,3,4,5,6};
				double[][] conditionExpects = new double[t + 1][];		
				IloLinearNumExpr sumdP = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) {
					Distribution[] distributions2 = new Distribution[1];
					distributions2[0] = new PoissonDist(sumLambdas[j][t]);
					
					PiecewiseComplementaryFirstOrderLossFunction pwcfolf = new PiecewiseComplementaryFirstOrderLossFunction(distributions2, seed);
					conditionExpects[j] = pwcfolf.getConditionalExpectations(probs, nbSamples);
					//System.out.println("max error is: " + pwcfolf.getMaxApproximationError(probs, nbSamples));
					sumdP.addTerm(P[j][t], sumLambdas[j][t]);
				}		
								
				IloLinearNumExpr sumPrdP = cplex.linearNumExpr();
				IloNumExpr xsumdP = cplex.linearNumExpr();
				for (int i = 0; i < segmentNum; i++) {
					double slope = 0; double interval = 0;	
					
					for (int k = 0; k <= i; k++) {
						slope += probs[k];
					}
					xsumdP= cplex.prod(cplex.sum(x[t], sumdP), slope);
				
					for (int j = 0; j <= t; j++) {
						for (int k = 0; k <= i; k++) {
							interval -= probs[k] * conditionExpects[j][k];
						}
						sumPrdP.addTerm(interval, P[j][t]);
					}
					cplex.addGe(xplus[t], cplex.sum(xsumdP, sumPrdP));
				}
				
				//inventory flow: d_t = S_t - I_t
				cplex.addEq(distributions[t].getMean(), cplex.diff(S[t], x[t]));
				
				//relation of inventory: x_t = x_t^+ - x_t^-
				cplex.addEq(x[t], cplex.diff(xplus[t], xminus[t]));
				
				if (t == 0) {
					tPurchaseCost = cplex.prod(v[t], cplex.diff(S[t], iniInventory));
					tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
					tTotalCost = cplex.sum(tHoldCost, cplex.sum(tTotalOrderCost, overheadCost));
									
					// cash flow: R_{t} = R_{0} + p_t(S_t - x^+_t) - hx^+_t - v(S_t - x^+_{t-1}) - Kz_t
					cplex.addEq(cplex.diff(R[t], iniCash), cplex.diff(tRevenue, tTotalCost));
					// cash constraint: R_{0} >= Kz_t + v(S_t - x^+_{t-1}) 
					cplex.addLe(cplex.sum(overheadCost, tTotalOrderCost), iniCash);					
				
					// z_t = 1 if S_1 - x^+_0 > 0
					// S_1 - x^+_0 >= 0
					// S_1 - x^+_0 <= z_t*M
					cplex.addGe(cplex.diff(S[t], iniInventory), 0);
					cplex.addLe(cplex.diff(S[t], iniInventory), cplex.prod(z[t], M));	
				}
				else {
					tPurchaseCost = cplex.prod(v[t], cplex.diff(S[t], xplus[t - 1]));
					tTotalOrderCost = cplex.sum(tPurchaseCost, tFixCost);
					tTotalCost = cplex.sum(tHoldCost, cplex.sum(tTotalOrderCost, overheadCost));
					
					// cash flow: R_{t} = R_{0} + p_t(S_t - x^+_t) - hx^+_t - v(S_t - x^+_{t-1}) - Kz_t
					cplex.addEq(cplex.diff(R[t], R[t - 1]), cplex.diff(tRevenue, tTotalCost)); 
					// cash constraint: R_{0} >= Kz_t + v(S_t - x^+_{t-1}) 
					cplex.addLe(cplex.sum(overheadCost, tTotalOrderCost), R[t - 1]); 

					// z_t = 1 if S_t - x^+_{t-1} > 0
					// S_t - x^+_{t-1} >= 0
					// S_t - x^+_{t-1} <= z_t*M
					cplex.addGe(cplex.diff(S[t], xplus[t-1]), 0); 
					cplex.addLe(cplex.diff(S[t], xplus[t-1]), cplex.prod(z[t], 10000)); 	
				}
				
				// x^+_t = max{S_t - d_t, 0}
				// x^+_t <= delta_t * M
				// S_t - d_t <= x^+_t + (1-delta_t)*M
				// S_t - d_t >= x^+_t - (1-delta_t)*M	
				cplex.addLe(xplus[t], cplex.prod(delta[t], 10000)); // s_t-I_{t-1} <= x_t*10000
				cplex.addLe(cplex.diff(S[t], distributions[t].getMean()), cplex.sum(xplus[t], cplex.prod(10000, cplex.diff(1, delta[t])))); 
				cplex.addGe(cplex.diff(S[t],distributions[t].getMean()), cplex.diff(xplus[t], cplex.prod(10000, cplex.diff(1, delta[t])))); 	
			}
						
			
			if (cplex.solve()) {				
				cplex.output().println("Solution status = " + cplex.getStatus());
				cplex.output().println("Solution value = " + cplex.getObjValue());
				System.out.println(cplex.getStatus());
				System.out.println("Solution value = " + cplex.getObjValue());
				varz = cplex.getValues(z);
				varS = cplex.getValues(S);
				varxplus = cplex.getValues(xplus);
				varx = cplex.getValues(x);
				varR = cplex.getValues(R);
				double[][] varP = new double[T][T];
				for (int i = 0; i < T; i++)
					for (int j = 0; j < T; j++)
						varP[i][j] = cplex.getValue(P[i][j]);
				System.out.println("z = ");
				System.out.println(Arrays.toString(varz));
				System.out.println("order-up-to level = ");
				System.out.println(Arrays.toString(varS));
				System.out.println("xplus = ");
				System.out.println(Arrays.toString(varxplus));
				System.out.println("x = ");
				System.out.println(Arrays.toString(varx));
				System.out.println("R = ");
				System.out.println(Arrays.toString(varR));
				System.out.println("P = ");
				System.out.println(Arrays.deepToString(varP));				
			}
			cplex.end();
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}

		return heuristicFindsCS(varz, varxplus, varR);
	}
	
	
	/**
	 * @param x
	 * @param varI
	 * @param varB
	 * @return find s, C(x), S based on the results of the MIP model
	 * @date: Nov 2, 2020, 11:07:07 AM 
	 */
	public double [][] heuristicFindsCS(double[] varx, double[] varI, double[] varB){
		// find s, C, S
		int T = distributions.length;
		double M = 10000;
		double[][] sCS = new double[T][3];

		// default values
		for (int t = 0; t < T; t++) {
			sCS[t][2] = distributions[T - 1].inverseF((price - variCost) / (holdingCost + price - salvageValue));
			sCS[t][1] = fixOrderCost;
		}

		double S = sCS[T - 1][2];
		double s = sCS[T - 1][1];
		for (int j = (int) S; j >= 0; j--) {
			if (Ly(j, T - 1, distributions[T - 1]) < Ly(S, T - 1, distributions[T - 1]) - fixOrderCost) {
				sCS[T - 1][0] = j + 1;
				break;
			}
		}
		sCS[T - 1][1] = 0; // C default value is 0
		// ascertain C for the last period
		for (int j = (int) S; j >= 0; j--) {
			int jj = 0;
			for (jj = j + 1; jj <= (int) S; jj++) {
				if (Ly(jj, T - 1, distributions[T - 1]) > fixOrderCost + Ly(j, T - 1, distributions[T - 1])) {
					sCS[T - 1][1] = fixOrderCost + variCost * (jj - 1 - j); // C for x = 0 at last period
					cacheC1Values.put(new State(T, j), sCS[T - 1][1]);
					break;
				}
			}
			if (Ly(S, T - 1, distributions[T - 1]) < fixOrderCost) { // choose a large value for C1, since expected
																		// profit is too small
				sCS[T - 1][1] = M;
				cacheC1Values.put(new State(T, j), sCS[T - 1][1]);
			}
		}

		for (int t = 0; t < T - 1; t++) {
			S = distributions[t].inverseF((price - variCost) / (holdingCost + price));
			if (Arrays.stream(varx).sum() < 0.1) // default values when x are all zeros
				break;
			int orderLastTo = nextIndex(varx, t);
			double demandSum = IntStream.range(t, orderLastTo + 1).mapToObj(i -> distributions[i].getMean()).reduce(0.0,
					(x, y) -> x.doubleValue() + y.doubleValue());
			double sigmaSqureSum = IntStream.range(t, orderLastTo + 1)
					.mapToObj(i -> distributions[i].getStandardDeviation())
					.reduce(0.0, (x, y) -> Math.pow(x.doubleValue(), 2) + Math.pow(y.doubleValue(), 2));
			double sigmaSum = Math.sqrt(sigmaSqureSum);
			Distribution distribution;
			if (distributions[0] instanceof ContinuousDistribution) // normal distribution
				distribution = new NormalDist(demandSum, sigmaSum);
			else
				distribution = new PoissonDist(demandSum);
			S = distribution.inverseF((price - variCost) / (holdingCost + price));

			// try a different S
			int orderLastTo2 = T - 1;
			double demandSum2 = IntStream.range(t, orderLastTo2 + 1).mapToObj(i -> distributions[i].getMean())
					.reduce(0.0, (x, y) -> x.doubleValue() + y.doubleValue());
			Distribution distribution2 = new PoissonDist(demandSum2);
			double S2 = distribution2.inverseF((price - variCost) / (price - salvageValue));

			double maxQ = 0;
			if (t == 0)
				maxQ = Math.max(0, (iniCash - fixOrderCost) / variCost);
			else
				maxQ = Math.max(0, (varB[t - 1] - fixOrderCost) / variCost);
			double cashS = t == 0 ? iniInventory + maxQ : varI[t - 1] + maxQ;
			S = Math.min(S, cashS);
			sCS[t][2] = S;
			//sCS[t][2] = S2; // can make optimality gap smaller when inventory holding cost is zero

			// ascertain s
			for (int i = 0; i <= (int) S; i++) {
				if (Ly(S, t, distribution) - Ly(i, t, distribution) < fixOrderCost + 0.1) {
					s = i;
					sCS[t][0] = s;
					break;
				}
				if (i == (int) S)
					sCS[t][0] = S;
			}
			if (sCS[0][0] == 0 && iniInventory == 0) {
				sCS[0][0] = 1;
			}

			// ascertain C
			for (int j = 0; j < (int) s; j++) {
				int jj = 0;
				S = distributions[t].inverseF((price - variCost) / (holdingCost + price)); // one period S
				for (jj = j + 1; jj <= (int) S; jj++) {
					if (Ly(jj, t, distributions[t]) > fixOrderCost + Ly(j, t, distributions[t])) {
						sCS[t][1] = fixOrderCost + variCost * (jj - 1 - j);
						cacheC1Values.put(new State(t + 1, j), sCS[t][1]);
						break;
					}
				}
				// System.out.println(Ly(onePeriodS, t, distributions[t]));
				if (Ly(S, t, distributions[t]) < fixOrderCost
						|| Ly(S, t, distributions[t]) - Ly(j, t, distributions[t]) < fixOrderCost // sometimes it is
																									// better to not use
																									// this condition
				) { // choose a large value for C, since expected profit is too small
					sCS[t][1] = fixOrderCost * 20;
					cacheC1Values.put(new State(t + 1, j), sCS[t][1]);
				}
			}

		}
		//sCS[1][0] = 1;
		// sCS[2][2] = 35; sCS[3][2] = 68;
		System.out.println("(s, C, S) from MIP are: " + Arrays.deepToString(sCS));
		return sCS;
	}
	
}
