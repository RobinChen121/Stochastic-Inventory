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
* @date: 2018年11月8日 下午9:44:23  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  this is a java class to implement the methods of Xiang and Rossi (2018) in ejor,
*                to obtain values of s and S by Joint mip model; when 4 periods, running time is 1s.
*                when test 8 periods, problem size limits exceed in CPLEX. It seems Ilog can tackle larger
*                size problems than cplex in java, but rather slow, more than 1 hour not get a solution.
*                
*/

public class JointMILP {
	double[] meanDemand; 
	double[] sigma;
	double fixOrderCost;
	double variCost;
	double holdingCost;
	double penaltyCost;
	double M;
	double[][] conSigma;
	int partionNum;
	BoundCriteria boundCriteria;
	boolean outputResults;
	public enum BoundCriteria{
		LOWBOUND,
		UPBOUND
	}
	
	public JointMILP(double[] meanDemand, double[] sigma, Double fixOrderCost, double variCost, double holdingCost,
			double penaltyCost, int partionNum, BoundCriteria boundCriteria) {
		this.meanDemand = meanDemand;
		this.sigma = sigma;
		this.fixOrderCost = fixOrderCost;
		this.variCost = variCost;
		this.holdingCost = holdingCost;
		this.penaltyCost = penaltyCost;
		this.partionNum = partionNum;
		this.M = 100000;
		this.boundCriteria = boundCriteria;
	}
	
	/*******************************************************************
	 * solve mip model by cplex, piecewise approximation
	 * @note: MipsS class must be initialized before invoking this method
	 */
	public double[] solveCPlex(int startPeriod) {
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
			break;
		default:
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
			IloIntVar[] xS = cplex.boolVarArray(T);  // whether ordering in period t
			IloNumVar[][] PS = new IloNumVar[T][T];
			for (int i = 0; i < PS.length; i++)
				for (int j = 0; j < PS.length; j++)
					PS[i][j] = cplex.boolVar();			
			IloNumVar[] IS = cplex.numVarArray(T, -Double.MAX_VALUE, Double.MAX_VALUE);
			IloNumVar[] IplusS = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // positive inventory
			IloNumVar[] IminusS = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // minus inventory
			
			IloIntVar[] xs = cplex.boolVarArray(T);  // whether ordering in period t
			IloNumVar[][] Ps = new IloNumVar[T][T];
			for (int i = 0; i < Ps.length; i++)
				for (int j = 0; j < Ps.length; j++)
					Ps[i][j] = cplex.boolVar();
			IloNumVar[] Is = cplex.numVarArray(T, -Double.MAX_VALUE, Double.MAX_VALUE);
			IloNumVar[] Ipluss = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // positive inventory
			IloNumVar[] Iminuss = cplex.numVarArray(T, 0.0, Double.MAX_VALUE); // minus inventory
			IloNumVar I0S = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE);
			IloNumVar I0s = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE);
			for (int t = 0; t < T; t++) {
				cplex.add(PS[t]);
				cplex.add(Ps[t]);
			}

			// objective function
			IloLinearNumExpr setupCostsS = cplex.linearNumExpr();
			IloLinearNumExpr holdCostsS = cplex.linearNumExpr();
			IloLinearNumExpr penaCostsS = cplex.linearNumExpr();
			IloNumExpr variCostsS = cplex.numExpr();
			setupCostsS.addTerms(xS, S);
			holdCostsS.addTerms(IplusS, h);
			penaCostsS.addTerms(IminusS, pai);
			variCostsS = cplex.prod(v[0], cplex.diff(IS[T - 1], I0S));
			IloNumExpr costsC = cplex.sum(setupCostsS, variCostsS, holdCostsS, penaCostsS);
			
			IloLinearNumExpr setupCostss = cplex.linearNumExpr();
			IloLinearNumExpr holdCostss = cplex.linearNumExpr();
			IloLinearNumExpr penaCostss = cplex.linearNumExpr();
			IloNumExpr variCostss = cplex.numExpr();
			setupCostss.addTerms(xs, S);
			holdCostss.addTerms(Ipluss, h);
			penaCostss.addTerms(Iminuss, pai);			
			variCostss = cplex.prod(v[0], cplex.diff(Is[T - 1], I0s));
			IloNumExpr costsG = cplex.sum(setupCostss, variCostss, holdCostss, penaCostss);
			
			IloLinearNumExpr objectiveCostS = cplex.linearNumExpr();
			IloLinearNumExpr objectiveCosts = cplex.linearNumExpr();
			for (int t = 0; t < T; t++) {
				objectiveCostS.addTerm(S[t], xS[t]);
				objectiveCostS.addTerm(h[t], IplusS[t]);
				objectiveCostS.addTerm(pai[t], IminusS[t]);
				if (t > 0) {
					objectiveCosts.addTerm(S[t], xs[t]);
					objectiveCosts.addTerm(h[t], Ipluss[t]);
					objectiveCosts.addTerm(pai[t], Iminuss[t]);
				}
			}
			
			cplex.addMinimize(cplex.sum(objectiveCostS, objectiveCosts));
			
			// constraints
			cplex.addEq(costsG, costsC);
			cplex.addGe(I0S, I0s);
			cplex.addEq(I0S, cplex.sum(newMeanDemand[0], IS[0]));
			//cplex.addEq(I0s, 15);
			
			// relationship between x_t and Q_t (I_t + d_t - I_{t-1} <= M*x_t)
			// Q_t >= 0
			for (int t = 0; t < T; t++) {
				if (t == 0) {
					cplex.addLe(cplex.sum(IS[t], cplex.diff(newMeanDemand[t], I0S)), cplex.prod(xS[t], M));
					cplex.addGe(cplex.sum(IS[t], newMeanDemand[t]), I0S);
					cplex.addLe(cplex.sum(Is[t], cplex.diff(newMeanDemand[t], I0s)), cplex.prod(xs[t], M));
					cplex.addGe(cplex.sum(Is[t], newMeanDemand[t]), I0s);
				}
				else {
					cplex.addLe(cplex.sum(cplex.sum(IS[t], newMeanDemand[t]), cplex.negative(IS[t - 1])), cplex.prod(xS[t], M));
					cplex.addGe(cplex.sum(IS[t], newMeanDemand[t]), IS[t - 1]);
					cplex.addLe(cplex.sum(cplex.sum(Is[t], newMeanDemand[t]), cplex.negative(Is[t - 1])), cplex.prod(xs[t], M));
					cplex.addGe(cplex.sum(Is[t], newMeanDemand[t]), Is[t - 1]);
				}				
			}
			
			// sum Pjt == 1
			IloLinearNumExpr sumPjt1;
			IloLinearNumExpr sumxjt1;
			for (int t = 0; t < T; t++) {	
				sumPjt1 = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) 
					sumPjt1.addTerm(1, PS[j][t]);
				cplex.addEq(sumPjt1, 1); // upper triangle
				for (int j = t + 1; j < T; j++)
					cplex.addEq(PS[j][t], 0);  // other Pjt = 0, or else cannot output values
			}
			
			IloLinearNumExpr sumPjt2;
			IloLinearNumExpr sumxjt2;
			for (int t = 0; t < T; t++) {	
				sumPjt2 = cplex.linearNumExpr();
				for (int j = 0; j <= t; j++) 
					sumPjt2.addTerm(1, PS[j][t]);
				cplex.addEq(sumPjt2, 1); // upper triangle
				for (int j = t + 1; j < T; j++)
					cplex.addEq(PS[j][t], 0);  // other Pjt = 0, or else cannot output values
			}

			// Pjt >= x_j - sum_{j+1}^{t}x_k
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					sumxjt1 = cplex.linearNumExpr();
					for (int k = j + 1; k <= t; k++)
						sumxjt1.addTerm(xS[k], 1);
					cplex.addGe(PS[j][t], cplex.diff(xS[j], sumxjt1));
				}
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					sumxjt2 = cplex.linearNumExpr();
					for (int k = j + 1; k <= t; k++)
						sumxjt2.addTerm(xS[k], 1);
					cplex.addGe(PS[j][t], cplex.diff(xS[j], sumxjt2));
				}
			
			// Pjt >= x_j - sum_{j+1}^{t}x_k
			// sum_{1}{t}x_k == 0 => P[0][t] == 1, this constraints are important for the extra piecewise constraints
			IloLinearNumExpr sumxjtS;
			IloLinearNumExpr sumxjts;
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					sumxjtS = cplex.linearNumExpr();
					sumxjts = cplex.linearNumExpr();
					for (int k = j + 1; k <= t; k++) {
						sumxjtS.addTerm(xS[k], 1);
						sumxjts.addTerm(xs[k], 1);
					}
					cplex.addGe(PS[j][t], cplex.diff(xS[j], sumxjtS));
					cplex.addGe(Ps[j][t], cplex.diff(xs[j], sumxjts));
					for (int k = 0; k <= t; k++) {
						sumxjtS.addTerm(xS[k], 1);
						sumxjts.addTerm(xs[k], 1);
					}
					cplex.addGe(cplex.prod(M, sumxjtS), cplex.prod(M, cplex.diff(1, PS[0][t])));
					cplex.addGe(cplex.prod(M, sumxjts), cplex.prod(M, cplex.diff(1, Ps[0][t])));
				}
			
			// special constraints for Cx and Gy
			cplex.addEq(xS[0], 1);
			cplex.addEq(xs[0], 0);
			
			// piecewise constraints
			IloNumExpr Ipk;
			IloLinearNumExpr PSigma;
			IloNumExpr pmeanPSigma;
			for (int t = 0; t < T; t++) {				
				for (int i = 0; i < partionNum; i++) {
					PSigma = cplex.linearNumExpr();
					double pik = Arrays.stream(prob).limit(i + 1).sum();
					Ipk = cplex.prod(IS[t], pik);
					
					double pmean = 0;
					for (int k = 0; k <= i; k++)
						pmean += prob[k] * means[k];
					
					for (int k = 0; k <= t; k++)
						PSigma.addTerm(PS[k][t], conSigma[k][t]);
									
					// upper bound					
					pmeanPSigma = cplex.prod(pmean, PSigma);
					IloNumExpr IpkMinuspmeanPSigma = cplex.diff(Ipk, pmeanPSigma);
					
					switch (boundCriteria) {
					case UPBOUND:
						// Iplus
						cplex.addGe(IplusS[t], cplex.sum(IpkMinuspmeanPSigma, cplex.prod(error, PSigma)));
						cplex.addGe(IplusS[t], cplex.prod(error, PSigma));
						
						// Iminus
						cplex.addGe(cplex.sum(IminusS[t], IS[t]), cplex.sum(IpkMinuspmeanPSigma, cplex.prod(error, PSigma)));
						cplex.addGe(cplex.sum(IminusS[t], IS[t]), cplex.prod(error, PSigma));
						break;
					
					case LOWBOUND:
						// Iplus
						cplex.addGe(IplusS[t], IpkMinuspmeanPSigma);
						cplex.addGe(IplusS[t], 0); // not necessary
						
						// Iminus
						cplex.addGe(cplex.sum(IminusS[t], IS[t]), IpkMinuspmeanPSigma);
						cplex.addGe(cplex.sum(IminusS[t], IS[t]), 0);
						break;					
					default:
						break;
					}
				}
			}
			
			Ipk = cplex.linearNumExpr();
			pmeanPSigma = cplex.linearNumExpr();
			for (int t = 0; t < T; t++) {				
				for (int i = 0; i < partionNum; i++) {
					PSigma = cplex.linearNumExpr();
					double pik = Arrays.stream(prob).limit(i + 1).sum();
					Ipk = cplex.prod(Is[t], pik);
					
					double pmean = 0;
					for (int k = 0; k <= i; k++)
						pmean += prob[k] * means[k];
					
					for (int k = 0; k <= t; k++)
						PSigma.addTerm(Ps[k][t], conSigma[k][t]);
									
					// upper bound					
					pmeanPSigma = cplex.prod(pmean, PSigma);
					IloNumExpr IpkMinuspmeanPSigma = cplex.diff(Ipk, pmeanPSigma);
					
					switch (boundCriteria) {
					case UPBOUND:
						// Iplus
						cplex.addGe(Ipluss[t], cplex.sum(IpkMinuspmeanPSigma, cplex.prod(error, PSigma)));
						cplex.addGe(Ipluss[t], cplex.prod(error, PSigma));
						
						// Iminus
						cplex.addGe(cplex.sum(Iminuss[t], Is[t]), cplex.sum(IpkMinuspmeanPSigma, cplex.prod(error, PSigma)));
						cplex.addGe(cplex.sum(Iminuss[t], Is[t]), cplex.prod(error, PSigma));
						break;
					
					case LOWBOUND:
						// Iplus
						cplex.addGe(Ipluss[t], IpkMinuspmeanPSigma);
						cplex.addGe(Ipluss[t], 0); // not necessary
						
						// Iminus
						cplex.addGe(cplex.sum(Iminuss[t], Is[t]), IpkMinuspmeanPSigma);
						cplex.addGe(cplex.sum(Iminuss[t], Is[t]), 0);
						break;					
					default:
						break;
					}
				}
			}
			
			// add another piecewise constraints, make results more robust
			IloNumExpr HMinusPiecewise;
			IloNumExpr BMinusPiecewise;
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					// Iplus
					double[] slopes = new double[partionNum + 1];
					double[] breakPointXCoor = means;
					slopes[0] = 0;
					for (int k = 1; k <= partionNum; k++)
						slopes[k] = slopes[k - 1] + prob[k - 1];
					double fa = 0;
					// require partionNum to be even, or else fa = 1/2 (prob[k]means[k] + prob[k+1]means[k+1]), k=(parNum-1)/2
					for(int k = 0; k < partionNum/2; k++)
						fa -= prob[k]*means[k];
					if(boundCriteria == BoundCriteria.UPBOUND)
						fa = error + fa;
					HMinusPiecewise = cplex.diff(cplex.prod(IplusS[t], 1/conSigma[j][t]), cplex.piecewiseLinear(cplex.prod(IS[t], 1/conSigma[j][t]), breakPointXCoor, slopes, 0, fa));
					cplex.addLe(HMinusPiecewise, cplex.prod(M, cplex.diff(1, PS[j][t])));
					cplex.addGe(HMinusPiecewise, cplex.prod(-M, cplex.diff(1, PS[j][t])));

					// Iminus
					slopes[0] = -1;
					for (int k = 1; k <= partionNum; k++)
						slopes[k] = slopes[k - 1] + prob[k - 1];
					BMinusPiecewise = cplex.diff(cplex.prod(IminusS[t], 1/conSigma[j][t]), cplex.piecewiseLinear(cplex.prod(IS[t], 1/conSigma[j][t]), breakPointXCoor, slopes, 0, fa));
					cplex.addLe(BMinusPiecewise, cplex.prod(M, cplex.diff(1, PS[j][t])));
					cplex.addGe(BMinusPiecewise, cplex.prod(-M, cplex.diff(1, PS[j][t])));
				}
			
			for (int t = 0; t < T; t++)
				for (int j = 0; j <= t; j++) {
					// Iplus
					double[] slopes = new double[partionNum + 1];
					double[] breakPointXCoor = means;
					slopes[0] = 0;
					for (int k = 1; k <= partionNum; k++)
						slopes[k] = slopes[k - 1] + prob[k - 1];
					double fa = 0;
					// require partionNum to be even, or else fa = 1/2 (prob[k]means[k] + prob[k+1]means[k+1]), k=(parNum-1)/2
					for(int k = 0; k < partionNum/2; k++)
						fa -= prob[k]*means[k];
					if(boundCriteria == BoundCriteria.UPBOUND)
						fa = error + fa;
					HMinusPiecewise = cplex.diff(cplex.prod(Ipluss[t], 1/conSigma[j][t]), cplex.piecewiseLinear(cplex.prod(Is[t], 1/conSigma[j][t]), breakPointXCoor, slopes, 0, fa));
					cplex.addLe(HMinusPiecewise, cplex.prod(M, cplex.diff(1, Ps[j][t])));
					cplex.addGe(HMinusPiecewise, cplex.prod(-M, cplex.diff(1, Ps[j][t])));

					// Iminus
					slopes[0] = -1;
					for (int k = 1; k <= partionNum; k++)
						slopes[k] = slopes[k - 1] + prob[k - 1];
					BMinusPiecewise = cplex.diff(cplex.prod(Iminuss[t], 1/conSigma[j][t]), cplex.piecewiseLinear(cplex.prod(Is[t], 1/conSigma[j][t]), breakPointXCoor, slopes, 0, fa));
					cplex.addLe(BMinusPiecewise, cplex.prod(M, cplex.diff(1, Ps[j][t])));
					cplex.addGe(BMinusPiecewise, cplex.prod(-M, cplex.diff(1, Ps[j][t])));
				}
			
			if (cplex.solve()) {				
//				double[] varx1 = cplex.getValues(xS);
//				double[] varI1 = cplex.getValues(IS);
//				double[] varIplus1 = cplex.getValues(IplusS);
//				double[] varIminus1 = cplex.getValues(IminusS);
//				double[] varx2 = cplex.getValues(xs);
//				double[] varI2 = cplex.getValues(Is);
//				double[] varIplus2 = cplex.getValues(Ipluss);
//				double[] varIminus2 = cplex.getValues(Iminuss);
				double[][] varP1 = new double[T][T];
				double[][] varP2 = new double[T][T];
				for (int i = 0; i < T; i++)
					for (int j = 0; j < T; j++) {
						varP1[i][j] = cplex.getValue(PS[i][j]);
						varP2[i][j] = cplex.getValue(Ps[i][j]);
					}
				
				double[] result = new double[2];
				System.out.println("Solution value = " + cplex.getObjValue());
				System.out.println("Solution status = " + cplex.getStatus());				
				if (outputResults = true) {
//					System.out.println("Solution status = " + cplex.getStatus());
//					System.out.println("s = " + cplex.getValue(I0s));
//					System.out.println("S = " + cplex.getValue(I0S));
//					System.out.println("costG is " + cplex.getValue(costsG));
//					System.out.println("costC is " + cplex.getValue(costsC));
//					System.out.println("xS = ");
//					System.out.println(Arrays.toString(varx1));
//					System.out.println("xs = ");
//					System.out.println(Arrays.toString(varx2));
//					System.out.println("IS = ");
//					System.out.println(Arrays.toString(varI1));
//					System.out.println("Is = ");
//					System.out.println(Arrays.toString(varI2));
//					String bound = boundCriteria == BoundCriteria.LOWBOUND ? "lower bound" : "upper bound";
//					System.out.println("IplusS " + bound + " = ");
//					System.out.println(Arrays.toString(varIplus1));
//					System.out.println("Ipluss " + bound + " = ");
//					System.out.println(Arrays.toString(varIplus2));
//					System.out.println("IminusS " + bound + "  = ");
//					System.out.println(Arrays.toString(varIminus1));
//					System.out.println("Iminuss " + bound + "  = ");
//					System.out.println(Arrays.toString(varIminus2));
//					System.out.println("PS = ");
//					System.out.println(Arrays.deepToString(varP1));	
//					System.out.println("Ps = ");
//					System.out.println(Arrays.deepToString(varP2));
					result[0] = cplex.getValue(I0s);
					result[1] = cplex.getValue(I0S);
					cplex.end();
				}
				return result;
			}			
			
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		
		return null;
	}
	
	public double[][] getsS() {
		int T = meanDemand.length;
		double[][] sS = new double[T][2];
		for (int t = 0; t < T; t++) {
			double[] sSt = solveCPlex(t);
			sS[t][0] = sSt[0];
			sS[t][1] = sSt[1];
		}
		return sS;
	}

	public static void main(String[] args) {
		double[] meanDemand = {10, 10,  10,	10,	10,	10,	10,	10};
		double[] sigma = Arrays.stream(meanDemand).map(i -> 0.25*i).toArray();
		double fixOrderCost = 100;
		double variCost = 0;
		double holdingCost = 1;
		double penaltyCost = 10;		
		int partionNum = 10;
		BoundCriteria boundCriteria = BoundCriteria.UPBOUND;
		
		long currTime = System.currentTimeMillis();
		JointMILP mipsS = new JointMILP(meanDemand, sigma, fixOrderCost, variCost, holdingCost, penaltyCost, partionNum, boundCriteria);
		double[][] sS = mipsS.getsS();
		System.out.println(Arrays.deepToString(sS));
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		/*******************************************************************
		 * Simulating Joint MILP
		 */
		int sampleNum = 10000;
		Simulation simulate = new Simulation(meanDemand, sigma, 0, fixOrderCost, variCost, holdingCost, penaltyCost);
		simulate.simulatesS(sS, sampleNum);

	}

}


