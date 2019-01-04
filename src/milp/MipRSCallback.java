package milp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.function.IntToDoubleFunction;
import java.util.stream.IntStream;

import org.apache.commons.math3.analysis.function.Abs;

import cern.jet.random.Normal;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLPMatrix;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloModeler;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.CplexStatus;
import sun.awt.www.content.audio.x_aiff;
import umontreal.ssj.functions.MathFunction;
import umontreal.ssj.functions.MathFunctionUtil;
import umontreal.ssj.probdist.NormalDist;

/**
* @author Zhen Chen
* @date: 2018年12月1日 上午9:51:52  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  This is a class to test call back method of Tuck et al. (2018) by the mip Rs approximation.
*                
* @note: this class need cplex.jar    
*        12 periods will exceed cplex default size.
*/

public class MipRSCallback {
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
	boolean outputResults;
	double[] cumSumDemand;
	
	public MipRSCallback(double[] meanDemand, double[] sigma, double iniInventory, Double fixOrderCost, double variCost, double holdingCost,
				double penaltyCost, int partionNum, boolean outputResults) {
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
		this.cumSumDemand = new double[T];
		this.cumSumDemand[0] = meanDemand[0];
		for (int t = 1; t < T; t++)
			this.cumSumDemand[t] = meanDemand[t] + cumSumDemand[t - 1];
		this.partionNum = partionNum;
		this.M = 100000;
		this.outputResults = outputResults;
	}
	
	// callback setttings
	public static class Callback extends IloCplex.UserCutCallback{
		IloRange[] cuts;
		
	    Callback(IloRange[] cuts) {
	    	this.cuts = cuts;
	    }
	    
	    public void main() throws IloException {
	    	int num = cuts.length;
	         for (int i = 0; i < num; ++i) {
	            IloRange thecut = cuts[i];
	            if ( thecut != null ) {
	            	add(thecut, IloCplex.CutManagement.UseCutForce); // add cuts 
	            	//cuts[i] = null;
	            }
	         }
	    }
	}
	
	public static IloRange[] makeCuts(IloCplex cplex, IloLPMatrix lp) 
												throws IloException{
		ArrayList<IloRange> cuts = new ArrayList<IloRange>();
		IloNumVar[] vars = lp.getNumVars();

		return cuts.toArray(new IloRange[cuts.size()]);
	}
	
	
	
	/*******************************************************************
	 * compute the expected left inventory for an ordering cycle from startPeriod to endPeriod,
	 * given order-up-to level S at startPeriod
	 */	
	double computeExpectIminus(double S, int startPeriod, int endPeriod) {
		double sigma = conSigma[startPeriod - 1][endPeriod - 1];
		double mu = IntStream.range(startPeriod - 1, endPeriod).mapToDouble(i -> meanDemand[i]).sum();
		MathFunction func = new MathFunction() {
			
			@Override
			public double evaluate(double x) {				
				return ((S - mu)/sigma - x) * NormalDist.density01(x);
			}
		};

		return sigma * MathFunctionUtil.integral(func, -50, (S - mu)/sigma) - S + mu; // I^-
		//return sigma * MathFunctionUtil.integral(func, -50, (S - mu)/sigma); // I^+
	}
	
	
	   
	/*******************************************************************
	 * solve mip model by cplex, piecewise approximation
	 * @note: MipRS class must be initialized before invoking this method
	 */
	public double solveCPlex() {
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
					for (int t = i; t <= j; t++) {
						Dxij.addTerm(-h[i] * cumSumDemand[t], x[i][j]);						
						holdCosts.addTerm(q[i][j], h[i]);
						//Dxij.addTerm(pai[i] * cumSumDemand[t], x[i][j]);
						//holdCosts.addTerm(q[i][j], -pai[i]);
						penaCosts.addTerm(h[i] + pai[i], H[i][j][t]);
					}
				}
			cplex.addMinimize(cplex.sum(setupCosts, holdCosts, Dxij, penaCosts));

			// constraints
			IloLPMatrix lp = cplex.addLPMatrix();
			
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
			lp.addRow(cplex.range(0, cplex.diff(sumX1j, 1), 0));
			lp.addRow(cplex.range(0, cplex.diff(sumXiT, 1), 0));
			for (int t = 0; t < T - 1; t++) {
				sumXit = cplex.linearNumExpr();
				sumXtj = cplex.linearNumExpr();
				for (int i = 0; i <= t; i++)
					sumXit.addTerm(1, x[i][t]);
				for (int j = t + 1; j < T; j++)
					sumXtj.addTerm(1, x[t + 1][j]);
				lp.addRow(cplex.range(0, cplex.diff(sumXit, sumXtj), 0));
			}

			// q_{i,j} <= Mx_{i,j}
			IloLinearNumExpr sumqit;
			IloLinearNumExpr sumqtj;
			for (int i = 0; i < T; i++)
				for (int j = 0; j < T; j++) 
					lp.addRow(cplex.range(-Double.MAX_VALUE, cplex.diff(q[i][j], cplex.prod(M, x[i][j])), 0));

			// sum_{i=0}^t q_{i, t} <= sum_{j=t+1}^T q_{t, j}			
			for (int t = 0; t < T - 1; t++) {
				sumqit = cplex.linearNumExpr();
				sumqtj = cplex.linearNumExpr();
				for (int i = 0; i <= t; i++)
					sumqit.addTerm(1, q[i][t]); 
				for (int j = t + 1; j < T; j++)
					sumqtj.addTerm(1, q[t + 1][j]);
				lp.addRow(cplex.range(-Double.MAX_VALUE, cplex.diff(sumqit, sumqtj), 0));
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
							
							IloNumExpr expr1 = cplex.sum(H[i][j][t], expectI);
							//IloNumExpr expr2 = cplex.diff(cplex.prod(expectI, pik), cplex.prod(x[i][j], conSigma[i][t] * pmean));
							//lp.addRow(cplex.range(0, cplex.diff(expr1, expr2), Double.MAX_VALUE));
							lp.addRow(cplex.range(0, expr1, Double.MAX_VALUE));							
						}
					}											
				}			
			
			cplex.setParam(IloCplex.Param.MIP.Strategy.Search, IloCplex.MIPSearch.Traditional);
			int iniConNum = lp.getNrows();
			System.out.println("number of constriants are:" + iniConNum);
			// maybe the slope and intercept error
			if (cplex.solve()) {	
				System.out.println("Solution value = " + cplex.getObjValue());
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
				
				System.out.println("z = ");
				System.out.println(Arrays.toString(z));
				System.out.println("Ordering quantities = ");
				System.out.println(Arrays.toString(quantity));
				System.out.println("I = ");
				System.out.println(Arrays.toString(I));
				System.out.println("x = ");
				System.out.println(Arrays.deepToString(varX));
				System.out.println("q = ");
				System.out.println(Arrays.deepToString(varQ));
				System.out.println("H = ");
				System.out.println(Arrays.deepToString(varH));
				
				double eps = 0.01;
				boolean addLine = true;
				ArrayList<IloRange> cuts = new ArrayList<>();
				while (addLine) {
					for (int i = 0; i < T; i++)
						for (int j = i; j < T; j++)
							for (int t = i; t <= j; t++) {
								if (cplex.getValue(x[i][j]) == 1) {
									double upToLevel = i > 0 ? cplex.getValue(q[i][j]) - cumSumDemand[i - 1] : cplex.getValue(q[i][j]);
									double Iminus = computeExpectIminus(upToLevel, i + 1, j + 1);								
									addLine = false;
									System.out.println("Iminus is " + Iminus);
									System.out.println("H[i][j][t] is " + cplex.getValue(H[i][j][t]));
									double valueH = cplex.getValue(H[i][j][t]);
									if (Iminus - valueH > eps) {
										double mu = IntStream.range(i, t + 1).mapToDouble(m -> meanDemand[m]).sum();
										double sigma = conSigma[i][t];
										double slope = NormalDist.cdf(mu, sigma, upToLevel) - 1;
										double intercept = Iminus - slope * upToLevel;
										IloNumExpr var = i > 0 ? cplex.diff(q[i][j], cplex.diff(x[i][j], cumSumDemand[i - 1]))
												: q[i][j];
										IloNumExpr tanLine = cplex.sum(intercept, cplex.prod(var, slope));
										IloRange cut = cplex.range(0, cplex.diff(H[i][j][t], tanLine), Double.MAX_VALUE);
										cuts.add(cut);
										addLine = true;									
									}
								}
							}
					IloRange[] cutsRange = cuts.toArray(new IloRange[cuts.size()]);
					lp.addRows(cutsRange);
					System.out.println("number of constriants are:" + lp.getNrows());
					cplex.solve();
					System.out.println("Solution status = " + cplex.getStatus());
					System.out.println("Solution value = " + cplex.getObjValue());
					System.out.println("number of constriants are:" + lp.getNrows());
					z = new double[T];
					quantity = new double[T];
					I = new double[T];
					lastQ = 0;
					for (int i = 0; i < T; i++)
						for (int j = 0; j < T; j++) {
							if (cplex.getValue(x[i][j]) == 1) {
								z[i] = 1;									
								if (i == 0) {
									quantity[i] = cplex.getValue(q[i][j]);
									lastQ = cplex.getValue(q[i][j]);
								}
								else
									quantity[i] = cplex.getValue(q[i][j]) - lastQ;
							}
						}
					I[0] = quantity[0] + iniInventory - meanDemand[0];
					for (int i = 1; i < T; i++)
						I[i] = quantity[i] + I[i - 1] - meanDemand[i];
					System.out.println("z = ");
					System.out.println(Arrays.toString(z));
					System.out.println("Ordering quantities = ");
					System.out.println(Arrays.toString(quantity));
					System.out.println("I = ");
					System.out.println(Arrays.toString(I));
					System.out.println("**************************************");
					System.out.println();
					lp.removeRows(iniConNum, cuts.size());
				}				
				
				if (outputResults == true) {					
//					System.out.println("x = ");
//					System.out.println(Arrays.deepToString(varX));
//					System.out.println("q = ");
//					System.out.println(Arrays.deepToString(varQ));
//					System.out.println("H = ");
//					System.out.println(Arrays.deepToString(varH));	
					
				}
				double finalOptValue = cplex.getObjValue();
				cplex.end();
				return finalOptValue;
			
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
		int partionNum = 10;
		boolean outputResults = true;
		long currTime = System.currentTimeMillis();
		MipRSCallback mipRS = new MipRSCallback(meanDemand, sigma, iniInventory, fixOrderCost, variCost, holdingCost, penaltyCost, partionNum, outputResults);
		mipRS.solveCPlex();
		double time = (System.currentTimeMillis() - currTime) / 1000;
		System.out.println("running time is " + time + "s");
		
		
		/*******************************************************************
		 * draw approximate picture for Gy 
		 */
//		int minInventorys = 0;
//		int maxInventorys = 200; // for drawing pictures
//		int xLength = maxInventorys - minInventorys + 1;
//		double[][] yG = new double[xLength][2];
//		int index = 0;
//		for (int  i = minInventorys; i <= maxInventorys; i++) {
//			iniInventory = i;
//			mipRS = new MipRS(meanDemand, sigma, iniInventory, fixOrderCost, variCost, holdingCost, penaltyCost, partionNum, boundCriteria, ComputeGyCx.COMPUTG, false);
//			yG[index][0] = i;
//			yG[index][1] = mipRS.solveCPlex();
//			index++;
//		}
//		Drawing.drawSimpleG(yG);
	}
}


