package milp;

import java.util.Arrays;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

/**
* @author Zhen Chen
* @date: 2018年11月14日 上午10:00:03  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  this is a class to implement the joint milp model in Xiang et al. (2018) on EJOR.
*/

public class JointMILP {
	double[] meanDemand; 
	double[] sigma;
	double iniInventory;	
	double fixOrderCost;
	double variCost;
	double holdingCost;
	double penaltyCost;
	double M;
	double[][] conSigma;
	int partionNum;
	BoundCriteria boundCriteria;
	public enum BoundCriteria{
		LOWBOUND,
		UPBOUND
	}
	
	public JointMILP(double[] meanDemand, double[] sigma, double iniInventory, Double fixOrderCost, double variCost, double holdingCost,
				double penaltyCost, int partionNum, BoundCriteria boundCriteria) {
		this.meanDemand = meanDemand;
		this.sigma = sigma;
		this.iniInventory = iniInventory;
		this.fixOrderCost = fixOrderCost;
		this.variCost = variCost;
		this.holdingCost = holdingCost;
		this.penaltyCost = penaltyCost;
		this.partionNum = partionNum;
		this.M = 100000;
		this.boundCriteria = boundCriteria;
	}
	
	/*******************************************************************
	 * solve mip model by cplex, get values of s by approximating C(x) and its relationship
	 * with G(y)
	 */
	public double[] solveCGets(double[] SG, int startPeriod) {
		// piecewise approximation values
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
		default:
			partionNum = 4;
			prob = new double[] {0.187555, 0.312445, 0.312445, 0.187555};
			means = new double[] {-1.43535, -0.415223, 0.415223, 1.43535};
			error = 0.0339052;
			break;
		}
		
		int T = meanDemand.length - startPeriod;
		double[] newMeanDemand = new double[T];
		double[] newSigma = new double[T];
		for (int i = 0; i < T; i++) {
			newMeanDemand[i] = meanDemand[i + startPeriod];
			newSigma[i] = sigma[i + startPeriod];
		}		
		
		conSigma = new double[T][T];
		for (int i = 0; i < T; i++)
			for (int j = 0; j < T; j++) {
				double sigmaPow = 0;
				for (int k = i; k <= j; k++) {
					sigmaPow += Math.pow(newSigma[k], 2);
				}
				conSigma[i][j] = Math.sqrt(sigmaPow);
			}
		
		try {
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); // no cplex logging information
			
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
			IloIntVar[] x = cplex.boolVarArray(T);  // whether ordering in period t
			IloNumVar[][] P = new IloNumVar[T][T];
			for (int i = 0; i < P.length; i++)
				for (int j = 0; j < P.length; j++)
					P[i][j] = cplex.boolVar();
			IloNumVar[] I = cplex.numVarArray(T, -Double.MAX_VALUE, Double.MAX_VALUE);
			IloNumVar[] Iplus = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // positive inventory
			IloNumVar[] Iminus = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // minus inventory
			//double I0 = iniInventory;
			IloNumVar I0 = cplex.numVar(0, Double.MAX_VALUE);
			
			// objective function
			IloLinearNumExpr setupCosts = cplex.linearNumExpr();
			IloLinearNumExpr holdCosts = cplex.linearNumExpr();
			IloLinearNumExpr penaCosts = cplex.linearNumExpr();
			IloNumExpr variCosts = cplex.numExpr();
			
			setupCosts.addTerms(x, S);
			holdCosts.addTerms(h, Iplus);
			penaCosts.addTerms(pai, Iminus);
			variCosts = cplex.prod(v[0], cplex.diff(I[T - 1], I0));
			IloNumExpr totalCosts = cplex.sum(setupCosts, variCosts, holdCosts, penaCosts);
			cplex.addMinimize(totalCosts);
			
			
			// constraints
			cplex.addEq(totalCosts, SG[1] + fixOrderCost);
			cplex.addLe(I0, SG[0]);
			
			//cplex.addEq(I0, 14);
			// relationship between x_t and Q_t (I_t + d_t - I_{t-1} <= M*x_t)
			// Q_t >= 0
			for (int t = 0; t < T; t++) {
				if (t == 0) {
					cplex.addLe(cplex.sum(I[t], cplex.diff(newMeanDemand[t], I0)), cplex.prod(x[t], M));
					cplex.addGe(cplex.sum(I[t], newMeanDemand[t]), I0);
				}
				else {
					cplex.addLe(cplex.sum(cplex.sum(I[t], newMeanDemand[t]), cplex.negative(I[t - 1])), cplex.prod(x[t], M));
					cplex.addGe(cplex.sum(I[t], newMeanDemand[t]), I[t - 1]);
				}				
			}
			
			// sum Pjt == 1
			IloLinearNumExpr sumPjt;
			IloLinearNumExpr sumxjt;
			for (int t = 0; t < T; t++) {	
				sumPjt = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) 
					sumPjt.addTerm(1, P[j][t]);
				cplex.addEq(sumPjt, 1); // upper triangle
				for (int j = t + 1; j < T; j++)
					cplex.addEq(P[j][t], 0);  // other Pjt = 0, or else cannot output values
			}
			
			// Pjt >= x_j - sum_{j+1}^{t}x_k
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					sumxjt = cplex.linearNumExpr();
					for (int k = j + 1; k <= t; k++)
						sumxjt.addTerm(x[k], 1);
					cplex.addGe(P[j][t], cplex.diff(x[j], sumxjt));
				}				
			
			//  piecewise constraints
			IloNumExpr Ipk;
			IloLinearNumExpr PSigma;
			IloNumExpr pmeanPSigma;
			for (int t = 0; t < T; t++) {				
				for (int i = 0; i < partionNum; i++) {
					PSigma = cplex.linearNumExpr();
					double pik = Arrays.stream(prob).limit(i + 1).sum();
					Ipk = cplex.prod(I[t], pik);
					
					double pmean = 0;
					for (int k = 0; k <= i; k++)
						pmean += prob[k] * means[k];
					
					for (int k = 0; k <= t; k++)
						PSigma.addTerm(P[k][t], conSigma[k][t]);
									
					// upper bound					
					pmeanPSigma = cplex.prod(pmean, PSigma);
					IloNumExpr IpkMinuspmeanPSigma = cplex.diff(Ipk, pmeanPSigma);
					
					switch (boundCriteria) {
					case UPBOUND:
						// Iplus
						cplex.addGe(Iplus[t], cplex.sum(IpkMinuspmeanPSigma, cplex.prod(error, PSigma)));
						cplex.addGe(Iplus[t], cplex.prod(error, PSigma));
						
						// Iminus
						cplex.addGe(cplex.sum(Iminus[t], I[t]), cplex.sum(IpkMinuspmeanPSigma, cplex.prod(error, PSigma)));
						cplex.addGe(cplex.sum(Iminus[t], I[t]), cplex.prod(error, PSigma));
						break;
					
					case LOWBOUND:
						// Iplus
						cplex.addGe(Iplus[t], IpkMinuspmeanPSigma);
						cplex.addGe(Iplus[t], 0); // not necessary
						
						// Iminus
						cplex.addGe(cplex.sum(Iminus[t], I[t]), IpkMinuspmeanPSigma);
						cplex.addGe(cplex.sum(Iminus[t], I[t]), 0);
						break;					
					default:
						break;
					}
				}
			}
			
			if (cplex.solve()) {	
				double[][] varP = new double[T][T];
				for (int i = 0; i < T; i++)
					for (int j = 0; j < T; j++)
						varP[i][j] = cplex.getValue(P[i][j]);
				
				System.out.println("Solution value = " + cplex.getObjValue());
				System.out.println("Solution status = " + cplex.getStatus());
				System.out.println("s = ");
				System.out.println(cplex.getValue(I0));
				double[] result = {cplex.getValue(I0), cplex.getObjValue()};
				return result;
			}
			cplex.end();
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		return null;
	}
	
	/*******************************************************************
	 * solve mip model by cplex, get values of S by approximating G(y)
	 */
	public double[] solveGGetS(int startPeriod) {
		// piecewise approximation values
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
		default:
			partionNum = 4;
			prob = new double[] {0.187555, 0.312445, 0.312445, 0.187555};
			means = new double[] {-1.43535, -0.415223, 0.415223, 1.43535};
			error = 0.0339052;
			break;
		}
		
		try {
			IloCplex cplex = new IloCplex();
			cplex.setOut(null); // no cplex logging information
			
			int T = meanDemand.length - startPeriod;
			double[] newMeanDemand = new double[T];
			double[] newSigma = new double[T];
			for (int i = 0; i < T; i++) {
				newMeanDemand[i] = meanDemand[i + startPeriod];
				newSigma[i] = sigma[i + startPeriod];
			}		
			
			conSigma = new double[T][T];
			for (int i = 0; i < T; i++)
				for (int j = 0; j < T; j++) {
					double sigmaPow = 0;
					for (int k = i; k <= j; k++) {
						sigmaPow += Math.pow(newSigma[k], 2);
					}
					conSigma[i][j] = Math.sqrt(sigmaPow);
				}
			
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
			IloIntVar[] x = cplex.boolVarArray(T);  // whether ordering in period t
			IloNumVar[][] P = new IloNumVar[T][T];
			for (int i = 0; i < P.length; i++)
				for (int j = 0; j < P.length; j++)
					P[i][j] = cplex.boolVar();
			IloNumVar[] I = cplex.numVarArray(T, -Double.MAX_VALUE, Double.MAX_VALUE);
			IloNumVar[] Iplus = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // positive inventory
			IloNumVar[] Iminus = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // minus inventory
			//double I0 = iniInventory;
			IloNumVar I0 = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE);
			
			// objective function
			IloLinearNumExpr setupCosts = cplex.linearNumExpr();
			IloLinearNumExpr holdCosts = cplex.linearNumExpr();
			IloLinearNumExpr penaCosts = cplex.linearNumExpr();
			IloNumExpr variCosts = cplex.numExpr();
			
			setupCosts.addTerms(x, S);
			holdCosts.addTerms(h, Iplus);
			penaCosts.addTerms(pai, Iminus);
			variCosts = cplex.prod(v[0], cplex.diff(I[T - 1], I0));
			cplex.addMinimize(cplex.sum(setupCosts, variCosts, holdCosts, penaCosts));
			
			
			// constraints
			
			// relationship between x_t and Q_t (I_t + d_t - I_{t-1} <= M*x_t)
			// Q_t >= 0
			for (int t = 0; t < T; t++) {
				if (t == 0) {
					cplex.addLe(cplex.sum(I[t], cplex.diff(newMeanDemand[t], I0)), cplex.prod(x[t], M));
					cplex.addGe(cplex.sum(I[t], newMeanDemand[t]), I0);
				}
				else {
					cplex.addLe(cplex.sum(cplex.sum(I[t], newMeanDemand[t]), cplex.negative(I[t - 1])), cplex.prod(x[t], M));
					cplex.addGe(cplex.sum(I[t], newMeanDemand[t]), I[t - 1]);
				}				
			}
			
			// sum Pjt == 1
			IloLinearNumExpr sumPjt;
			IloLinearNumExpr sumxjt;
			for (int t = 0; t < T; t++) {	
				sumPjt = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) 
					sumPjt.addTerm(1, P[j][t]);
				cplex.addEq(sumPjt, 1); // upper triangle
				for (int j = t + 1; j < T; j++)
					cplex.addEq(P[j][t], 0);  // other Pjt = 0, or else cannot output values
			}
			
			// Pjt >= x_j - sum_{j+1}^{t}x_k
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					sumxjt = cplex.linearNumExpr();
					for (int k = j + 1; k <= t; k++)
						sumxjt.addTerm(x[k], 1);
					cplex.addGe(P[j][t], cplex.diff(x[j], sumxjt));
				}
			
			// for computing G(y)
			cplex.addEq(x[0], 0);		
			
			//  piecewise constraints
			IloNumExpr Ipk;
			IloLinearNumExpr PSigma;
			IloNumExpr pmeanPSigma;
			for (int t = 0; t < T; t++) {				
				for (int i = 0; i < partionNum; i++) {
					PSigma = cplex.linearNumExpr();
					double pik = Arrays.stream(prob).limit(i + 1).sum();
					Ipk = cplex.prod(I[t], pik);
					
					double pmean = 0;
					for (int k = 0; k <= i; k++)
						pmean += prob[k] * means[k];
					
					for (int k = 0; k <= t; k++)
						PSigma.addTerm(P[k][t], conSigma[k][t]);
									
					// upper bound					
					pmeanPSigma = cplex.prod(pmean, PSigma);
					IloNumExpr IpkMinuspmeanPSigma = cplex.diff(Ipk, pmeanPSigma);
					
					switch (boundCriteria) {
					case UPBOUND:
						// Iplus
						cplex.addGe(Iplus[t], cplex.sum(IpkMinuspmeanPSigma, cplex.prod(error, PSigma)));
						cplex.addGe(Iplus[t], cplex.prod(error, PSigma));
						
						// Iminus
						cplex.addGe(cplex.sum(Iminus[t], I[t]), cplex.sum(IpkMinuspmeanPSigma, cplex.prod(error, PSigma)));
						cplex.addGe(cplex.sum(Iminus[t], I[t]), cplex.prod(error, PSigma));
						break;
					
					case LOWBOUND:
						// Iplus
						cplex.addGe(Iplus[t], IpkMinuspmeanPSigma);
						cplex.addGe(Iplus[t], 0); // not necessary
						
						// Iminus
						cplex.addGe(cplex.sum(Iminus[t], I[t]), IpkMinuspmeanPSigma);
						cplex.addGe(cplex.sum(Iminus[t], I[t]), 0);
						break;					
					default:
						break;
					}
				}
			}
			
			if (cplex.solve()) {								
				System.out.println("Solution value = " + cplex.getObjValue());
				System.out.println("Solution status = " + cplex.getStatus());
				System.out.println("S = " + cplex.getValue(I0));
				double[] result = {cplex.getValue(I0), cplex.getObjValue()};
				return result;
			}
			cplex.end();
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		return null;
	}
	
	/*******************************************************************
	 * get values of s and S for each period by linear approximation
	 */
	public double[][] getAppsS() {
		int T = meanDemand.length;
		double[][] sS = new double[T][2];
		for (int t = 0; t < T; t++) {
			double[] SG = solveGGetS(t);
			double[] sC = solveCGets(SG, t);
			sS[t][0] = sC[0];
			sS[t][1] = SG[0];
		}
		return sS;
	}

	public static void main(String[] args) {
		double[] meanDemand = {20, 40, 60, 40};
		double[] sigma = Arrays.stream(meanDemand).map(i -> 0.25*i).toArray();
		double iniInventory = 0;	
		double fixOrderCost = 100;
		double variCost = 0;
		double holdingCost = 1;
		double penaltyCost = 10;		
		int partionNum = 10;
		BoundCriteria boundCriteria = BoundCriteria.LOWBOUND;
		
		long currTime = System.currentTimeMillis();
		JointMILP joint = new JointMILP(meanDemand, sigma, iniInventory, fixOrderCost, variCost, holdingCost, penaltyCost, partionNum, boundCriteria);
		double[][] sS = joint.getAppsS();
		System.out.println(Arrays.deepToString(sS));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		/*******************************************************************
		 * Simulating Joint MILP
		 */
		int sampleNum = 10000;
		Simulation simulate = new Simulation(meanDemand, sigma, iniInventory, fixOrderCost, variCost, holdingCost, penaltyCost);
		simulate.simulatesS(sS, sampleNum);

	}

}


