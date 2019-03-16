package milp;

import java.util.Arrays;
import java.util.function.IntFunction;
import java.util.stream.IntStream;


import ilog.concert.IloConstraint;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import umontreal.ssj.stochprocess.GeometricNormalInverseGaussianProcess;


/**
* @author Zhen Chen
* @date: 2018年11月28日 下午6:38:53  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  this is a class to implement the PM setting method of Tunc et al. (2018) to solve 
*                single-item stochastic lot sizing problem.
*                (1) without dynamic cut, cplex reach size limit even for 8 periods
*                (2) the reason of using H[i][j][t] may lie on the convenience of adopting call back
*/

public class MipRSPM {
	double[] meanDemand; 
	double[] sigma;
	double iniInventory;	
	double fixOrderCost;
	private double variCost;
	private double holdingCost;
	private double penaltyCost;
	double[] cumSumDemand;
	int partionNum;
	int T;
	int M;
	boolean outputResults;
	double[][] conSigma;
	
	public MipRSPM(double[] meanDemand, double[] sigma, double iniInventory, Double fixOrderCost, double variCost, double holdingCost,
			double penaltyCost, int partionNum, boolean outputResults) {
		this.meanDemand = meanDemand;
		this.sigma = sigma;
		this.iniInventory = iniInventory;
		this.fixOrderCost = fixOrderCost;
		this.variCost = variCost;
		this.holdingCost = holdingCost;
		this.penaltyCost = penaltyCost;
		this.T = meanDemand.length;
		this.M = 100000;	
		this.cumSumDemand = new double[T];
		this.cumSumDemand[0] = meanDemand[0];
		for (int t = 1; t < T; t++)
			this.cumSumDemand[t] = meanDemand[t] + cumSumDemand[t - 1];
		this.partionNum = partionNum;
		this.outputResults = outputResults;
		this.conSigma = new double[T][T];
		for (int i = 0; i < T; i++)
			for (int j = 0; j < T; j++) {
				double sigmaPow = 0;
				for (int k = i; k <= j; k++) {
					sigmaPow += Math.pow(sigma[k], 2);
				}
				this.conSigma[i][j] = Math.sqrt(sigmaPow);
			}
	}
	
	public double solveCallBack() {
		// piecewise approximation values
		// H \approx aS+b, a = \sum_{i=1}^w prob[i], b = \sum_{i=1}^w prob[i]*mean[i], multiply sigma[i][j] for general norm
		// plus error for upper bound
		double[] prob;
		double[] means;
		double error;
		switch (partionNum) {
		case 4:
			prob = new double[] {0.187555, 0.312445, 0.312445, 0.187555};
			means = new double[] {-1.43535, -0.415223, 0.415223, 1.43535};
			error = 0.0339052;
			break;

		case 10:
			prob = new double[] {0.04206108420763477, 0.0836356495308449, 0.11074334596058821, 0.1276821455299152, 0.13587777477101692, 0.13587777477101692, 0.1276821455299152, 0.11074334596058821, 0.0836356495308449, 0.04206108420763477};
			means = new double[] {-2.133986195498256, -1.3976822972668839, -0.918199946431143, -0.5265753462727588, -0.17199013069262026, 0.17199013069262026, 0.5265753462727588, 0.918199946431143, 1.3976822972668839, 2.133986195498256};
			error = 0.005885974956458359;
			break;
		default:
			prob = new double[] {0.187555, 0.312445, 0.312445, 0.187555};
			means = new double[] {-1.43535, -0.415223, 0.415223, 1.43535};
			error = 0.0339052;
			break;
		}
				
		try {
			IloCplex cplex = new IloCplex();
			cplex.getVersion();
			cplex.setOut(null); // no cplex logging information
			
			// parameter values in array
			double[] K = new double[T];
			double[] h = new double[T];
			double[] v = new double[T];
			double[] pai = new double[T];
			Arrays.fill(K, fixOrderCost);
			Arrays.fill(h, holdingCost);
			Arrays.fill(v, variCost);
			Arrays.fill(pai, penaltyCost);
			
			// decision variables			
			IloIntVar[][] x = new IloIntVar[T][T];   // whether [i, j) is a replenishment cycle
			IloNumVar[][] q = new IloNumVar[T][T];   // total ordering quantity sums to period i, if [i, j) is a replenishment cycle
			IloNumVar[][][] H = new IloNumVar[T][T][T]; // loss function value at period t of repre cycle [i, j)
			for (int i = 0; i < T; i++)
				for (int j = 0; j < T; j++) {					
					for (int t = 0; t < T; t++) {
						if (t >= i && t <= j)
							H[i][j][t] = cplex.numVar(0, Double.MAX_VALUE);
						else
							H[i][j][t] = cplex.numVar(0, 0);
					}
					if (j < i) {
						x[i][j] = cplex.intVar(0, 0);
						q[i][j] = cplex.numVar(0, 0);
					}
					else {
						x[i][j] = cplex.boolVar();
						q[i][j] = cplex.numVar(0, Double.MAX_VALUE);
					}
				}
			
			// objective function
			IloLinearNumExpr setupCosts = cplex.linearNumExpr();
			IloLinearNumExpr holdCosts = cplex.linearNumExpr();
			IloLinearNumExpr Dxij = cplex.linearNumExpr();
			IloLinearNumExpr penaCosts = cplex.linearNumExpr();
			for (int i = 0; i < T; i++)
				for (int j = 0; j < T; j++) {
					setupCosts.addTerm(x[i][j], K[i]);	
					
					for (int t = i; t <= j; t++) {
						//Dxij.addTerm(pai[i] * cumSumDemand[t], x[i][j]);
						//holdCosts.addTerm(q[i][j], -pai[i]);
						Dxij.addTerm(-h[i] * cumSumDemand[t], x[i][j]);
						holdCosts.addTerm(q[i][j], h[i]);
						penaCosts.addTerm(h[i] + pai[i], H[i][j][t]);
					}
				}
			cplex.addMinimize(cplex.sum(setupCosts, holdCosts, Dxij, penaCosts));
			
			// constraints
			// sum_{i=0}^T x_{1, i} = 1
			// sum_{i=0}^T x_{i, T} = 1;
			// sum_{i=0}^t x_{i, t} = sum_{j=t+1}^T x_{t, j}						
			IloLinearNumExpr sumX1j = cplex.linearNumExpr();
			IloLinearNumExpr sumXiT = cplex.linearNumExpr();
			IloLinearNumExpr sumXit;
			IloLinearNumExpr sumXtj;
			for (int i = 0; i < T; i++) {
				sumX1j.addTerm(1, x[0][i]);
				sumXiT.addTerm(1, x[i][T - 1]);
			}
			cplex.addEq(sumX1j, 1);
			cplex.addEq(sumXiT, 1);
			for (int t = 0; t < T - 1; t++) {
				sumXit = cplex.linearNumExpr();
				sumXtj = cplex.linearNumExpr();
				for (int i = 0; i <= t; i++)
					sumXit.addTerm(1, x[i][t]);
				for (int j = t + 1; j < T; j++)
					sumXtj.addTerm(1, x[t + 1][j]);
				cplex.addEq(sumXit, sumXtj);
			}
			
			// q_{i,j} <= Mx_{i,j}
			IloLinearNumExpr sumqit;
			IloLinearNumExpr sumqtj;
			for (int i = 0; i < T; i++)
				for (int j = 0; j < T; j++) 
					cplex.addLe(q[i][j], cplex.prod(M, x[i][j]));
			
			// sum_{i=0}^t q_{i, t} <= sum_{j=t+1}^T q_{t, j}			
			for (int t = 0; t < T - 1; t++) {
				sumqit = cplex.linearNumExpr();
				sumqtj = cplex.linearNumExpr();
				for (int i = 0; i <= t; i++)
					sumqit.addTerm(1, q[i][t]); 
				for (int j = t + 1; j < T; j++)
					sumqtj.addTerm(1, q[t + 1][j]);
				cplex.addLe(sumqit, sumqtj);
			}
			
			// piecewise constraints
			IloNumExpr expectI;
			for (int i = 0; i < T; i++)
				for (int j = i; j < T; j++) {		
					for (int t = i; t <= j; t++) {
						expectI = cplex.diff(q[i][j], cplex.prod(cumSumDemand[t], x[i][j]));						
						for (int k = 0; k < partionNum; k++) {
							double pik = Arrays.stream(prob).limit(k + 1).sum();							
							double pmean = 0;
							for (int m = 0; m <= k; m++)
								pmean += prob[m] * means[m];
							//cplex.addGe(H[i][j][t], cplex.diff(cplex.prod(expectI, pik), cplex.prod(x[i][j], conSigma[i][t] * pmean)));
							cplex.addGe(cplex.sum(H[i][j][t], expectI), cplex.diff(cplex.prod(expectI, pik), cplex.prod(x[i][j], conSigma[i][t] * pmean)));
						}
					}											
				}			
			 
			if (cplex.solve()) {				
				double[][][] varH = new double[T][T][T];
				double[][] varQ = new double[T][T];
				double[][] varX = new double[T][T];
				for (int i = 0; i < T; i++)
					for (int j = 0; j < T; j++) {
						varX[i][j] = cplex.getValue(x[i][j]);
						varQ[i][j] = cplex.getValue(q[i][j]);
						for (int k = i; k <= j; k++)
							varH[i][j][k] = cplex.getValue(H[i][j][k]);
					}
				double[] z = new double[T];
				double[] quantity = new double[T];
				double[] I = new double[T];
				double lastQ = 0;
				for (int i = 0; i < T; i++)
					for (int j = 0; j < T; j++) {
						if (varX[i][j] == 1) {
							z[i] = 1;								
							if (i == 0) {
								quantity[i] = varQ[i][j];
								lastQ = quantity[i];
							}
							else {
								quantity[i] = varQ[i][j] - lastQ;
								lastQ = quantity[i];
							}
						}
					}
				I[0] = quantity[0] + iniInventory - meanDemand[0];
				for (int i = 1; i < T; i++)
					I[i] = quantity[i] + I[i - 1] - meanDemand[i];
				
				System.out.println("Solution value = " + cplex.getObjValue());
				if (outputResults == true) {
					System.out.println("Solution status = " + cplex.getStatus());
					System.out.println("number of constriants are:" + cplex.getNcols()); // it's not less than number of variables
					System.out.println("z = ");
					System.out.println(Arrays.toString(z));
					System.out.println("Ordering quantities = ");
					System.out.println(Arrays.toString(quantity));
					System.out.println("I = ");
					System.out.println(Arrays.toString(I));
//					System.out.println("part of holding costs is :");
//					System.out.println(cplex.getValue(holdCosts));
//					System.out.println("Another part of holding costs is :");
//					System.out.println(cplex.getValue(Dxij));
//					System.out.println("merged penalty costs is :");
//					System.out.println(cplex.getValue(penaCosts));
					System.out.println("x = ");
					System.out.println(Arrays.deepToString(varX));
					System.out.println("q = ");
					System.out.println(Arrays.deepToString(varQ));
					System.out.println("H = ");
					System.out.println(Arrays.deepToString(varH));
	
				}
				double finalOptValue = cplex.getObjValue();
				cplex.end();
				return finalOptValue;
			}
			else {

			}
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		
		return 0;
	} 

	public static void main(String[] args) {
		double[] meanDemand = {50,50,50,50,50,50};
		double[] sigma = Arrays.stream(meanDemand).map(i -> 0.25*i).toArray();
		double iniInventory = 0;	
		double fixOrderCost = 100;
		double variCost = 0;
		double holdingCost = 1;
		double penaltyCost = 10;
		int partionNum = 10;
		boolean outputResults = true;
		
		MipRSPM mipCallBack = new MipRSPM(meanDemand, sigma, iniInventory, fixOrderCost, variCost, holdingCost,
				penaltyCost, partionNum, outputResults);
		mipCallBack.solveCallBack();

	}

}


