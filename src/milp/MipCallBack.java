package milp;

import java.util.Arrays;
import java.util.function.IntFunction;
import java.util.stream.IntStream;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

/**
* @author Zhen Chen
* @date: 2018年11月28日 下午6:38:53  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  this is a class to implement the calling back method of Tunc et al. (2018) to solve 
*                single-item stochastic lot sizing problem.
*/

public class MipCallBack {
	double[] meanDemand; 
	double[] sigma;
	double iniInventory;	
	double fixOrderCost;
	double variCost;
	double holdingCost;
	double penaltyCost;
	double[] cumSumDemand;
	int partionNum;
	int T;
	int M;
	
	public MipCallBack(double[] meanDemand, double[] sigma, double iniInventory, Double fixOrderCost, double variCost, double holdingCost,
			double penaltyCost, int partionNum) {
		this.meanDemand = meanDemand;
		this.sigma = sigma;
		this.iniInventory = iniInventory;
		this.fixOrderCost = fixOrderCost;
		this.variCost = variCost;
		this.holdingCost = holdingCost;
		this.penaltyCost = penaltyCost;
		this.T = meanDemand.length;
		this.M = 100000;	
		this.cumSumDemand[0] = meanDemand[0];
		for (int t = 1; t < T; t++)
			this.cumSumDemand[t] = meanDemand[t] + cumSumDemand[t - 1];
		this.partionNum = partionNum;
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
					for (int t = i; t < j; t++) {
						Dxij.addTerm(-cumSumDemand[t], x[i][j]);
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
				for (int j = t; j < T - 1; j++)
					sumXtj.addTerm(1, x[t][j]);
				cplex.addEq(sumXit, sumXtj);
			}
			
			// q_{i,j} <= Mx_{i,j}
			IloLinearNumExpr sumqit;
			IloLinearNumExpr sumqtj;
			for (int i = 0; i < T; i++)
				for (int j = 0; j < T; j++) 
					cplex.addLe(q[i][j], cplex.prod(M, x[i][j]));
			
			// sum_{i=0}^t q_{i, t} = sum_{j=t+1}^T q_{t, j}			
			for (int t = 0; t < T - 1; t++) {
				sumqit = cplex.linearNumExpr();
				sumqtj = cplex.linearNumExpr();
			}
			
			// piecewise constraints
			for (int i = 0; i < T; i++)
				for (int j = i; j < T; j++) {
					for (int t = i; t <= j; t++) {
						                                                                                 
					}
						
				}
			
			
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		
		return 0;
	} 

	public static void main(String[] args) {
		double[] meanDemand = {20, 40, 60, 40};
		double[] sigma = Arrays.stream(meanDemand).map(i -> 0.25*i).toArray();
		double iniInventory = 0;	
		double fixOrderCost = 100;
		double variCost = 0;
		double holdingCost = 1;
		double penaltyCost = 10;
		
		

	}

}


