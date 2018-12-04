package milp;

import java.util.Arrays;
import java.util.function.IntFunction;
import java.util.function.IntToDoubleFunction;
import java.util.stream.IntStream;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBExpr;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;




/**
* @author Zhen Chen
* @date: 2018年12月3日 下午4:14:34  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  this is class to implement the piecewise method of Roberto Rossi et al. (2015) in Gurobi of java
* 
*                Gurobi does not support very well for java. 
*                There are much fewer classes and methods for java compared with that of cplex
*/

public class MipRSGurobi {
	double[] meanDemand; 
	double[] sigma;
	double iniInventory;	
	double fixOrderCost;
	double variCost;
	double holdingCost;
	double penaltyCost;
	int T;
	double M;
	double[][] conSigma;
	int partionNum;
	BoundCriteria boundCriteria;
	ComputeGyCx gyCx;
	boolean outputResults;
	public enum BoundCriteria{
		LOWBOUND,
		UPBOUND
	}
	
	public enum ComputeGyCx{
		COMPUTG,
		COMPUTC,
		NOTCOMPUT;
	}
	
	public MipRSGurobi(double[] meanDemand, double[] sigma, double iniInventory, Double fixOrderCost, double variCost, double holdingCost,
			double penaltyCost, int partionNum, BoundCriteria boundCriteria, ComputeGyCx gyCx, boolean outputResults) {
	this.meanDemand = meanDemand;
	this.sigma = sigma;
	this.iniInventory = iniInventory;
	this.fixOrderCost = fixOrderCost;
	this.variCost = variCost;
	this.holdingCost = holdingCost;
	this.penaltyCost = penaltyCost;
	this.T = meanDemand.length;
	conSigma = new double[T][T];
	for (int i = 0; i < T; i++)
		for (int j = 0; j < T; j++) {
			double sigmaPow = 0;
			for (int k = i; k <= j; k++) {
				sigmaPow += Math.pow(sigma[k], 2);
			}
			this.conSigma[i][j] = Math.sqrt(sigmaPow);
		}
	this.partionNum = partionNum;
	this.M = 100000;
	this.boundCriteria = boundCriteria;
	this.gyCx = gyCx;
	this.outputResults = outputResults;
}
	
	/*******************************************************************
	 * solve mip model by cplex, piecewise approximation
	 * @note: MipRS class must be initialized before invoking this method
	 */
	public double solveGurobi() {
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
			GRBEnv env = new GRBEnv();
			GRBModel model = new GRBModel(env);
			
			// parameter values in array
			double[] S = new double[T];
			double[] h = new double[T];
			double[] v = new double[T];
			double[] pai = new double[T];
			Arrays.fill(S, fixOrderCost);
			Arrays.fill(h, holdingCost);
			Arrays.fill(v, variCost);
			Arrays.fill(pai, penaltyCost);
			
			// decision variables	
			double I0 = iniInventory;
			GRBVar[] x = new GRBVar[T]; 
			GRBVar[] I = new GRBVar[T];
			GRBVar[] Iplus = new GRBVar[T];
			GRBVar[] Iminus = new GRBVar[T];
			GRBVar[][] P = new GRBVar[T][T];
			for (int t = 0; t < T; t++) {
				x[t] = model.addVar(0, 1, 0, GRB.BINARY, "x" + t);
				I[t] = model.addVar(-Double.MAX_VALUE, Double.MAX_VALUE, 0, GRB.CONTINUOUS, "I" + t);
				Iplus[t] = model.addVar(0, Double.MAX_VALUE, 0, GRB.CONTINUOUS, "Iplus" + t);
				Iminus[t] = model.addVar(0, Double.MAX_VALUE, 0, GRB.CONTINUOUS, "Iminus" + t);
				for (int j = 0; j < P.length; j++) {
					if (j <= t)
						P[t][j] = model.addVar(0, 1, 0, GRB.BINARY, "P" + t + j);
					else
						P[t][j] = model.addVar(0, 0, 0, GRB.BINARY, "P" + t + j);
				}
			}
			
			// objective function
			GRBLinExpr totalCosts = new GRBLinExpr();
			GRBLinExpr fixOrderCosts = new GRBLinExpr();
			GRBLinExpr variCosts = new GRBLinExpr();
			GRBLinExpr holdCosts = new GRBLinExpr();
			GRBLinExpr penaCosts = new GRBLinExpr();
			GRBLinExpr ItMinusI0 = new GRBLinExpr();
			
			fixOrderCosts.addTerms(S, x);
			holdCosts.addTerms(h, Iplus);
			penaCosts.addTerms(pai, Iminus); 
			ItMinusI0.addTerm(-I0, I[T - 1]);
			variCosts.multAdd(v[0], ItMinusI0);			
			totalCosts.add(holdCosts); totalCosts.add(penaCosts); totalCosts.add(variCosts); totalCosts.add(fixOrderCosts);
			
			// constraints
			
			// relationship between x_t and Q_t (I_t + d_t - I_{t-1} <= M*x_t)
			// Q_t >= 0
			GRBLinExpr expr1 = new GRBLinExpr();
			GRBLinExpr expr2 = new GRBLinExpr();
			GRBLinExpr expr3 = new GRBLinExpr();
			for (int t = 0; t < T; t++) {
				if (t == 0) {
					expr1.addTerm(meanDemand[t] - I0, I[t]);
					expr2.addTerm(M, x[t]);
					expr3.addTerm(meanDemand[t], I[t]);
					model.addConstr(expr1, GRB.LESS_EQUAL, expr2, "");
					model.addConstr(expr3, GRB.GREATER_EQUAL, I0, "");
				}
				else {
					
				}				
			}
			
			
			
		} catch (GRBException e) {
			 System.out.println("Error code: " + e.getErrorCode() + ". " +
                     e.getMessage());
		}
		return 0;
	}

	public static void main(String[] args) {
		double[] meanDemand = {20, 40, 60, 40, 20, 40, 60, 40};
		double[] sigma = Arrays.stream(meanDemand).map(i -> 0.25*i).toArray();
		double iniInventory = 0;	
		double fixOrderCost = 100;
		double variCost = 0;
		double holdingCost = 1;
		double penaltyCost = 10;		
		int partionNum = 10;

	}

}


