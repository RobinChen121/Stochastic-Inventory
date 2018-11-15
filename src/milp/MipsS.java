package milp;

import java.util.Arrays;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import milp.MipRS.BoundCriteria;
import sun.security.provider.JavaKeyStore.CaseExactJKS;


/**
* @author Zhen Chen
* @date: 2018年11月8日 下午9:44:23  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  this is a java class to implement the methods of Xiang and Rossi (2018) in ejor,
*                to obtain values of s and S by mip model;
*                the joint milp model seems couldn't obtain values of s and S directly
*/

public class MipsS {
	double[] meanDemand; 
	double[] sigma;
	double fixOrderCost;
	double variCost;
	double holdingCost;
	double penaltyCost;
	int T;
	double M;
	double[][] conSigma;
	int partionNum;
	BoundCriteria boundCriteria;
	boolean outputResults;
	public enum BoundCriteria{
		LOWBOUND,
		UPBOUND
	}
	
	public MipsS(double[] meanDemand, double[] sigma, Double fixOrderCost, double variCost, double holdingCost,
			double penaltyCost, int partionNum, BoundCriteria boundCriteria) {
		this.meanDemand = meanDemand;
		this.sigma = sigma;
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
	}
	
	/*******************************************************************
	 * solve mip model by cplex, piecewise approximation
	 * @note: MipsS class must be initialized before invoking this method
	 */
	public double solveCPlex() {
		// piecewise approximation values
		double[] prob;
		double[] means;
		double error;
		switch (partionNum) {
		case 4: // in fact, it is 5 segments
			prob = new double[] {0.187555, 0.312445, 0.312445, 0.187555};
			means = new double[] {-1.43535, -0.415223, 0.415223, 1.43535};
			error = 0.0339052;
			break;

		case 10: // in fact, it is 11 segments
			prob = new double[] {0.04206108420763477, 0.0836356495308449, 0.11074334596058821, 0.1276821455299152, 0.13587777477101692, 0.13587777477101692, 0.1276821455299152, 0.11074334596058821, 0.0836356495308449, 0.04206108420763477};
			means = new double[] {-2.133986195498256, -1.3976822972668839, -0.918199946431143, -0.5265753462727588, -0.17199013069262026, 0.17199013069262026, 0.5265753462727588, 0.918199946431143, 1.3976822972668839, 2.133986195498256};
			error = 0.005885974956458359;
			
		default:
			partionNum = 4;
			prob = new double[] {0.187555, 0.312445, 0.312445, 0.187555};
			means = new double[] {-1.43535, -0.415223, 0.415223, 1.43535};
			error = 0.0339052;
			break;
		}
		
		try {
			IloCplex cplex = new IloCplex();
			//cplex.setOut(null); // no cplex logging information
			
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
			IloIntVar[] x1 = cplex.boolVarArray(T);  // whether ordering in period t
			IloNumVar[][] P1 = new IloNumVar[T][T];
			for (int i = 0; i < P1.length; i++)
				for (int j = 0; j < P1.length; j++)
					P1[i][j] = cplex.boolVar();			
			IloNumVar[] I1 = cplex.numVarArray(T, -Double.MAX_VALUE, Double.MAX_VALUE);
			IloNumVar[] Iplus1 = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // positive inventory
			IloNumVar[] Iminus1 = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // minus inventory
			
			IloIntVar[] x2 = cplex.boolVarArray(T);  // whether ordering in period t
			IloNumVar[][] P2 = new IloNumVar[T][T];
			for (int i = 0; i < P2.length; i++)
				for (int j = 0; j < P2.length; j++)
					P2[i][j] = cplex.boolVar();
			IloNumVar[] I2 = cplex.numVarArray(T, -Double.MAX_VALUE, Double.MAX_VALUE);
			IloNumVar[] Iplus2 = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // positive inventory
			IloNumVar[] Iminus2 = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // minus inventory
			IloNumVar I01 = cplex.numVar(-100, 100);
			IloNumVar I02 = cplex.numVar(-100, 100);
			cplex.add(x1);
			cplex.add(x2);
			cplex.add(Iplus1);
			cplex.add(Iplus2);
			cplex.add(Iminus1);
			cplex.add(Iminus2);
			cplex.add(I01);
			cplex.add(I02);
			cplex.add(I1);
			cplex.add(I2);
			for (int t = 0; t < T; t++) {
				cplex.add(P1[t]);
				cplex.add(P2[t]);
			}

			// objective function
			IloLinearNumExpr setupCosts1 = cplex.linearNumExpr();
			IloLinearNumExpr holdCosts1 = cplex.linearNumExpr();
			IloLinearNumExpr penaCosts1 = cplex.linearNumExpr();
			IloNumExpr variCosts1 = cplex.numExpr();
			setupCosts1.addTerms(x1, S);
			holdCosts1.addTerms(Iplus1, h);
			penaCosts1.addTerms(Iminus1, pai);
			variCosts1 = cplex.prod(v[0], cplex.diff(I1[T - 1], I01));
			IloNumExpr costsC = cplex.sum(setupCosts1, variCosts1, holdCosts1, penaCosts1);
			
			IloLinearNumExpr setupCosts2 = cplex.linearNumExpr();
			IloLinearNumExpr holdCosts2 = cplex.linearNumExpr();
			IloLinearNumExpr penaCosts2 = cplex.linearNumExpr();
			IloNumExpr variCosts2 = cplex.numExpr();
			setupCosts2.addTerms(x2, S);
			holdCosts2.addTerms(Iplus2, h);
			penaCosts2.addTerms(Iminus2, pai);			
			variCosts2 = cplex.prod(v[0], cplex.diff(I2[T - 1], I02));
			IloNumExpr costsG = cplex.sum(setupCosts2, variCosts2, holdCosts2, penaCosts2);
			
			IloLinearNumExpr objectiveCost1 = cplex.linearNumExpr();
			IloLinearNumExpr objectiveCost2 = cplex.linearNumExpr();
			for (int t = 0; t < T; t++) {
				objectiveCost1.addTerm(S[t], x1[t]);
				objectiveCost1.addTerm(h[t], Iplus1[t]);
				objectiveCost1.addTerm(pai[t], Iplus1[t]);
				if (t > 0) {
					objectiveCost2.addTerm(S[t], x2[t]);
					objectiveCost2.addTerm(h[t], Iplus2[t]);
					objectiveCost2.addTerm(pai[t], Iplus2[t]);
				}
			}
			
			cplex.addMinimize(cplex.sum(objectiveCost1, objectiveCost2));
			
			// constraints
			cplex.addEq(costsG, costsC);
			cplex.addLe(I01, I02);
			cplex.addEq(I02, cplex.sum(I2[0], meanDemand[0]));
			
			// relationship between x_t and Q_t (I_t + d_t - I_{t-1} <= M*x_t)
			// Q_t >= 0
			for (int t = 0; t < T; t++) {
				if (t == 0) {
					cplex.addLe(cplex.sum(I1[t], cplex.diff(meanDemand[t], I01)), cplex.prod(x1[t], M));
					cplex.addGe(cplex.sum(I1[t], meanDemand[t]), I01);
					cplex.addLe(cplex.sum(I2[t], cplex.diff(meanDemand[t], I02)), cplex.prod(x2[t], M));
					cplex.addGe(cplex.sum(I2[t], meanDemand[t]), I02);
				}
				else {
					cplex.addLe(cplex.sum(cplex.sum(I1[t], meanDemand[t]), cplex.negative(I1[t - 1])), cplex.prod(x1[t], M));
					cplex.addGe(cplex.sum(I1[t], meanDemand[t]), I1[t - 1]);
					cplex.addLe(cplex.sum(cplex.sum(I2[t], meanDemand[t]), cplex.negative(I2[t - 1])), cplex.prod(x2[t], M));
					cplex.addGe(cplex.sum(I2[t], meanDemand[t]), I2[t - 1]);
				}				
			}
			
			// sum Pjt == 1
			IloLinearNumExpr sumPjt1;
			IloLinearNumExpr sumxjt1;
			for (int t = 0; t < T; t++) {	
				sumPjt1 = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) 
					sumPjt1.addTerm(1, P1[j][t]);
				cplex.addEq(sumPjt1, 1); // upper triangle
				for (int j = t + 1; j < T; j++)
					cplex.addEq(P1[j][t], 0);  // other Pjt = 0, or else cannot output values
			}
			
			IloLinearNumExpr sumPjt2;
			IloLinearNumExpr sumxjt2;
			for (int t = 0; t < T; t++) {	
				sumPjt2 = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) 
					sumPjt2.addTerm(1, P1[j][t]);
				cplex.addEq(sumPjt2, 1); // upper triangle
				for (int j = t + 1; j < T; j++)
					cplex.addEq(P1[j][t], 0);  // other Pjt = 0, or else cannot output values
			}

			// Pjt >= x_j - sum_{j+1}^{t}x_k
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					sumxjt1 = cplex.linearNumExpr();
					for (int k = j + 1; k <= t; k++)
						sumxjt1.addTerm(x1[k], 1);
					cplex.addGe(P1[j][t], cplex.diff(x1[j], sumxjt1));
				}
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					sumxjt2 = cplex.linearNumExpr();
					for (int k = j + 1; k <= t; k++)
						sumxjt2.addTerm(x1[k], 1);
					cplex.addGe(P1[j][t], cplex.diff(x1[j], sumxjt2));
				}
			
			// special constraints for Cx and Gy
			cplex.addEq(x1[0], 1);
			cplex.addEq(x2[0], 0);
			
		//  piecewise constraints
			IloNumExpr Ipk1;
			IloLinearNumExpr PSigma1;
			IloNumExpr pmeanPSigma1;
			IloNumExpr Ipk2;
			IloLinearNumExpr PSigma2;
			IloNumExpr pmeanPSigma2;
			for (int t = 0; t < T; t++) {				
				for (int i = 0; i < partionNum; i++) {
					PSigma1 = cplex.linearNumExpr();
					double pik1 = Arrays.stream(prob).limit(i + 1).sum();
					Ipk1 = cplex.prod(I1[t], pik1);
					double pmean1 = 0;
					for (int k = 0; k <= i; k++)
						pmean1 += prob[k] * means[k];
					for (int k = 0; k <= t; k++)
						PSigma1.addTerm(P1[k][t], conSigma[k][t]);
					
					PSigma2 = cplex.linearNumExpr();
					double pik2 = Arrays.stream(prob).limit(i + 1).sum();
					Ipk2 = cplex.prod(I2[t], pik2);
					double pmean2 = 0;
					for (int k = 0; k <= i; k++)
						pmean2 += prob[k] * means[k];
					for (int k = 0; k <= t; k++)
						PSigma2.addTerm(P2[k][t], conSigma[k][t]);

					// upper bound					
					pmeanPSigma1 = cplex.prod(pmean1, PSigma1);
					IloNumExpr IpkMinuspmeanPSigma1 = cplex.diff(Ipk1, pmeanPSigma1);
					
					pmeanPSigma2 = cplex.prod(pmean2, PSigma2);
					IloNumExpr IpkMinuspmeanPSigma2 = cplex.diff(Ipk2, pmeanPSigma2);

					switch (boundCriteria) {
					case UPBOUND:
						// Iplus
						cplex.addGe(Iplus1[t], cplex.sum(IpkMinuspmeanPSigma1, cplex.prod(error, PSigma1)));
						cplex.addGe(Iplus1[t], cplex.prod(error, PSigma1));
						// Iminus
						cplex.addGe(cplex.sum(Iminus1[t], I1[t]), cplex.sum(IpkMinuspmeanPSigma1, cplex.prod(error, PSigma1)));
						cplex.addGe(cplex.sum(Iminus1[t], I1[t]), cplex.prod(error, PSigma1));
						
						// Iplus
						cplex.addGe(Iplus2[t], cplex.sum(IpkMinuspmeanPSigma2, cplex.prod(error, PSigma2)));
						cplex.addGe(Iplus2[t], cplex.prod(error, PSigma2));
						// Iminus
						cplex.addGe(cplex.sum(Iminus2[t], I2[t]), cplex.sum(IpkMinuspmeanPSigma2, cplex.prod(error, PSigma2)));
						cplex.addGe(cplex.sum(Iminus2[t], I2[t]), cplex.prod(error, PSigma2));					
						break;
					case LOWBOUND:
						// Iplus
						cplex.addGe(Iplus1[t], IpkMinuspmeanPSigma1);
						cplex.addGe(Iplus1[t], 0); // not necessary
						// Iminus
						cplex.addGe(cplex.sum(Iminus1[t], I1[t]), IpkMinuspmeanPSigma1);
						cplex.addGe(cplex.sum(Iminus1[t], I1[t]), 0);
						
						// Iplus
						cplex.addGe(Iplus2[t], IpkMinuspmeanPSigma2);
						cplex.addGe(Iplus2[t], 0); // not necessary
						// Iminus
						cplex.addGe(cplex.sum(Iminus2[t], I2[t]), IpkMinuspmeanPSigma2);
						cplex.addGe(cplex.sum(Iminus2[t], I2[t]), 0);						
						break;					
					default:
						break;
					}
				}
			}
			
			if (cplex.solve()) {				
				double[] varx1 = cplex.getValues(x1);
				double[] varI1 = cplex.getValues(I1);
				double[] varIplus1 = cplex.getValues(Iplus1);
				double[] varIminus1 = cplex.getValues(Iminus1);
				double[][] varP1 = new double[T][T];
				double[] varx2 = cplex.getValues(x2);
				double[] varI2 = cplex.getValues(I2);
				double[] varIplus2 = cplex.getValues(Iplus2);
				double[] varIminus2 = cplex.getValues(Iminus2);
				double[][] varP2 = new double[T][T];
				for (int i = 0; i < T; i++)
					for (int j = 0; j < T; j++) {
						varP1[i][j] = cplex.getValue(P1[i][j]);
						varP2[i][j] = cplex.getValue(P2[i][j]);
					}
				
				System.out.println("Solution value = " + cplex.getObjValue());
				System.out.println("Solution status = " + cplex.getStatus());
				if (outputResults != true) {
					System.out.println("Solution status = " + cplex.getStatus());
					System.out.println("s = " + I02);
					System.out.println("S = " + I01);
					System.out.println("costG is " + cplex.getValue(costsG));
					System.out.println("costC is " + cplex.getValue(costsC));
					System.out.println("x1 = ");
					System.out.println(Arrays.toString(varx1));
					System.out.println("x2 = ");
					System.out.println(Arrays.toString(varx2));
					System.out.println("I1 = ");
					System.out.println(Arrays.toString(varI1));
					System.out.println("I2 = ");
					System.out.println(Arrays.toString(varI2));
					String bound = boundCriteria == BoundCriteria.LOWBOUND ? "lower bound" : "upper bound";
					System.out.println("Iplus1 " + bound + " = ");
					System.out.println(Arrays.toString(varIplus1));
					System.out.println("Iplus2 " + bound + " = ");
					System.out.println(Arrays.toString(varIplus2));
					System.out.println("Iminus1 " + bound + "  = ");
					System.out.println(Arrays.toString(varIminus1));
					System.out.println("Iminus2 " + bound + "  = ");
					System.out.println(Arrays.toString(varIminus2));
					System.out.println("P1 = ");
					System.out.println(Arrays.deepToString(varP1));	
					System.out.println("P2 = ");
					System.out.println(Arrays.deepToString(varP2));	
				}
				return cplex.getObjValue();
			}
			cplex.end();
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		
		return 0;
	}

	public static void main(String[] args) {
		double[] meanDemand = {40};
		double[] sigma = Arrays.stream(meanDemand).map(i -> 0.25*i).toArray();
		double fixOrderCost = 100;
		double variCost = 0;
		double holdingCost = 1;
		double penaltyCost = 10;		
		int partionNum = 10;
		BoundCriteria boundCriteria = BoundCriteria.UPBOUND;
		
		MipsS mipsS = new MipsS(meanDemand, sigma, fixOrderCost, variCost, holdingCost, penaltyCost, partionNum, boundCriteria);
		mipsS.solveCPlex();

	}

}


