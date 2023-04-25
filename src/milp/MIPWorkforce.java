package milp;


import static org.junit.jupiter.api.DynamicTest.stream;

import java.util.Arrays;

import org.apache.poi.hssf.record.LeftMarginRecord;
import org.apache.poi.poifs.crypt.DataSpaceMapUtils.IRMDSTransformInfo;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBQConstr;
import gurobi.GRBQuadExpr;
import gurobi.GRBVar;
import umontreal.ssj.probdist.BinomialDist;



public class MIPWorkforce {
	int iniStaffNum;
	double fixCost;
	double unitVariCost;
	double salary;
	double unitPenalty;		
	int[] minStaffNum;
	double[] turnoverRate;
	int T;
	
	public MIPWorkforce(int iniStaffNum, double fixCost, double unitVariCost, double salary, double unitPenalty,
			int[] minStaffNum, double[] turnoverRate) {
		this.iniStaffNum = iniStaffNum;
		this.fixCost = fixCost;
		this.unitVariCost = unitVariCost;
		this.unitPenalty = unitPenalty;
		this.salary = salary;
		this.minStaffNum = minStaffNum;		
		this.turnoverRate = turnoverRate;
		this.T = turnoverRate.length;
	}
	
	/**
	 * original formulation of the loss function for the workforce problem
	 * @param y
	 * @return
	 */
	public double lossFunction(int y, int w, double p) {
		BinomialDist dist = new BinomialDist(y, p);
		double value = 0;
		int imin = Math.max(y - w, 0); 
		for (int i = imin; i <= y; i++) {
			value += dist.prob(i) * (double) (i+w-y);
		}
		return value;
	}
	
	public double[][] piecewise(int segmentNum, int w, double p) {
		double[] slope = new double[segmentNum + 1];
		double[] intercept = new double[segmentNum + 1];
		double[] tanPointXcoe = new double[segmentNum + 1];
		double[] tanPointYcoe = new double[segmentNum + 1];
		double[] intsPointXcoe = new double[segmentNum + 2];
		double[] intsPointYcoe = new double[segmentNum + 2];
		double[] intsPointGap = new double[segmentNum + 2];
		double[][] result = new double[7][];
		
		int endX = w;
		for (int k = w + 1; k < w*10; k ++) {
			endX = k;
			if (Fy(k, w, p) > 0.9999) {
				endX = k;
				break;
			}
		}		
		
		slope[segmentNum] = 0;
		tanPointXcoe[segmentNum] = endX;
		tanPointYcoe[segmentNum] = 0;
		intercept[segmentNum] = 0;
		for (int i = 0; i < segmentNum; i++) {
			if (i == 0) {
				slope[i] = p - 1;
				tanPointXcoe[0] = w - 1;
				tanPointYcoe[0] = (w - 1) * p + 1;
				intercept[0] = w;
			}
			else {
				int a = (int)tanPointXcoe[i-1];
				tanPointXcoe[i] = a;
				slope[i] = slope[i-1];
				for (int j = a + 1; j <= endX; j++) {
					if (Fy(j, w, p) - Fy(a, w, p) > 1 /(double)segmentNum) {
						tanPointXcoe[i] = j;
						int b = (int)tanPointXcoe[i];
						tanPointYcoe[i] = lossFunction(b, w, p);
						slope[i] = -(1 - p)*(1 - Fy(b, w, p));
						intercept[i] = -slope[i] * tanPointXcoe[i] + tanPointYcoe[i];
						break;
					}
				}				
			}
		}
		intsPointXcoe[0] = 0;
		intsPointYcoe[0] = w * p;
		intsPointGap[0] = 0;
		intsPointXcoe[segmentNum + 1] = endX;
		intsPointYcoe[segmentNum + 1] = 0;
		intsPointGap[segmentNum + 1] = 0;
		for (int i = 0; i < segmentNum; i++) {
			intsPointXcoe[i+1] = tanPointYcoe[i+1] - tanPointYcoe[i] + slope[i]*tanPointXcoe[i] - slope[i+1]*tanPointXcoe[i+1];
			intsPointXcoe[i+1] = slope[i] ==  slope[i+1] ? tanPointXcoe[i] : intsPointXcoe[i+1] / (slope[i] -  slope[i+1]);
			intsPointYcoe[i+1] = slope[i]*(intsPointXcoe[i+1] - tanPointXcoe[i]) + tanPointYcoe[i];
			double y =lossFunction((int)intsPointXcoe[i+1], w, p);
			intsPointGap[i+1] = y - intsPointYcoe[i+1];
		}
		result[0] = slope;
		result[1] = intercept;
		result[2] = tanPointXcoe;
		result[3] = tanPointYcoe;
		result[4] = intsPointXcoe;
		result[5] = intsPointYcoe;
		result[6] = intsPointGap;
		return result;
	}
	
	/**
	 * @param y
	 * @return the cdf F_y(y-W)
	 */
	public double Fy(int y, int w, double p) {
		BinomialDist dist = new BinomialDist(y, p);
		return dist.cdf(y - w);
	}
	
	
	public double pieceApprox(int segmentNum) {		
		try {
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);		
			//env.set("logFile", "mip-chance.log");
			env.start();
			
			// Create empty model
			GRBModel model = new GRBModel(env);			
			model.set(GRB.IntParam.LogToConsole, 0); // disable console logging
			
			// Create variables
		    GRBVar[] y = new GRBVar[T];
		    GRBVar[] u = new GRBVar[T];
		    GRBVar[] x = new GRBVar[T];
		    GRBVar[] z = new GRBVar[T];
		    GRBVar[][] P = new GRBVar[T][T];
		    for (int t = 0; t < T; t++) {
		    	y[t] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "y");
		    	u[t] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "u");
		    	//delta[t] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "delta");
		    	x[t] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "x");
		    	z[t] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "z");
		    	for (int j = 0; j <= t; j++) {
		    		P[j][t] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "P");
		    	}
		    } 
		    
		   // objective function, set objective
		    GRBLinExpr obj = new GRBLinExpr();
		    for (int t = 0; t < T; t++) {
		    	obj.addTerm(fixCost, z[t]);
		    	if (t == 0) {
		    		obj.addTerm(unitVariCost, y[t]);
		    		obj.addConstant(-unitVariCost*iniStaffNum);
		    		
		    	}
		    	else {
		    		GRBLinExpr iExpr = new GRBLinExpr();
		    		iExpr.addTerm(1, y[t]);
		    		iExpr.addTerm(-1, x[t-1]);
		    		obj.multAdd(unitVariCost, iExpr);
				}
		    	obj.addTerm(unitPenalty, u[t]);
		    	obj.addTerm(salary, x[t]);
		    }
		    model.setObjective(obj, GRB.MINIMIZE);

		    
		    // constraints
		    int M = Integer.MAX_VALUE;
		    for (int t = 0; t < T; t++) {	
		    	GRBLinExpr left1 = new GRBLinExpr();	
		    	
		    	// y_t >= x_{t-1}
		    	// y_t - x_{t-1} <= z_t M
		    	if (t == 0) 
		    		left1.addConstant(-iniStaffNum);	
		    	else 
		    		left1.addTerm(-1, x[t-1]);	
		    	left1.addTerm(1, y[t]);
		    	model.addConstr(left1, GRB.GREATER_EQUAL, 0, null); 
		    	GRBLinExpr right1 = new GRBLinExpr();
		    	right1.addTerm(M, z[t]);
		    	model.addConstr(left1, GRB.LESS_EQUAL, right1, null);
		    	
		    	// sum_{j=1}^t P_{jt} == 1
		    	GRBLinExpr left2 = new GRBLinExpr();
		    	for (int j = 0; j <= t; j++)
		    		left2.addTerm(1, P[j][t]);
		    	model.addConstr(left2, GRB.EQUAL, 1, null);
		    	
//		    	model.addConstr(P[1][1], GRB.EQUAL, 1, null);
//		    	model.addConstr(y[1], GRB.EQUAL, 66, null);
//		    	model.addConstr(y[0], GRB.EQUAL, 66, null);
		    	
		    	// P_{jt} <= z_j - \sum_{k=j+1}^t z[k]
		    	for (int j = 0; j <= t; j++) {
		    		GRBLinExpr right2 = new GRBLinExpr();
		    		for (int k = j + 1; k <= t; k++) 
		    			right2.addTerm(-1, z[k]);
		    		right2.addTerm(1, z[j]);
		    		model.addConstr(P[j][t], GRB.LESS_EQUAL, right2, null);
		    	}
		    		    	
		    	// x_t >= y_j(1-p)^{t-j+1} - (1-P_{jt})M
		    	// x_t <= y_j(1-p)^{t-j+1} + (1-P_{jt})M
		    	for (int j = 0; j <= t; j++) {
		    		double p = 1;
		    		for (int k = j; k <= t; k++)
		    			p = p * (1 - turnoverRate[j]);
		    		GRBLinExpr right3 = new GRBLinExpr();
		    		right3.addTerm(p, y[j]);
		    		right3.addTerm(M, P[j][t]);
		    		right3.addConstant(-M);
		    		model.addConstr(x[t], GRB.GREATER_EQUAL, right3, null);
		    		GRBLinExpr right4 = new GRBLinExpr();
		    		right4.addTerm(p, y[j]);
		    		right4.addTerm(-M, P[j][t]);
		    		right4.addConstant(M);
		    		model.addConstr(x[t], GRB.LESS_EQUAL, right4, null);
		    	}
		    	
		    	// piecewise constraints
		    	// U_t >= \alpha y_j + \beta - (1 - P_{jt})M
		    	for (int j = 0; j <= t; j++) {
		    		double p = 1;
		    		for (int k = j; k <= t; k++)
		    			p = p * (1 - turnoverRate[j]);
		    		double[][] result = piecewise(segmentNum, minStaffNum[t], 1 - p);
		    		double[] slope = result[0];
		    		double[] intercept = result[1];
		    		double[] gap = result[6];
		    		double error = Arrays.stream(gap).max().getAsDouble();
		    		for (int m = 0; m < segmentNum; m++) {
		    			GRBLinExpr right5 = new GRBLinExpr();
		    			right5.addTerm(slope[m], y[j]);
		    			right5.addConstant(intercept[m]);
		    			
		    			// lower bound
		    			right5.addTerm(M, P[j][t]);
		    			right5.addConstant(-M);
		    			model.addConstr(u[t], GRB.GREATER_EQUAL, right5, null);
		    			
		    			// upper bound
//		    			right5.addTerm(M, P[j][t]);
//		    			right5.addConstant(-M);
//		    			right5.addConstant(error);		
//		    			model.addConstr(u[t], GRB.GREATER_EQUAL, right5, null);
		    		}
		    	}
		    }
		    
		    // Optimize model
			model.optimize();
			    
			// output results
			System.out.println(model.get(GRB.DoubleAttr.ObjVal));
			double[] yV = new double[T];
			double[] xV = new double[T];
			double[] uV = new double[T];
			double[] zV = new double[T];
			double PV;
			for (int t = 0; t < T; t++) {
				yV[t] = y[t].get(GRB.DoubleAttr.X);
				xV[t] = x[t].get(GRB.DoubleAttr.X);
				uV[t] = u[t].get(GRB.DoubleAttr.X);
				zV[t] = z[t].get(GRB.DoubleAttr.X);
			}
			PV = P[0][0].get(GRB.DoubleAttr.X);
			System.out.println("P is " + PV);
			System.out.println("z is " + Arrays.toString(zV));
			System.out.println("y is " + Arrays.toString(yV));
			System.out.println("x is " + Arrays.toString(xV));
			System.out.println("u is " + Arrays.toString(uV));
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
		}
		return 0;
	}
	

	public int[] getZ() {
		int[] Z = new int[T];
		try {	
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);		
			//env.set("logFile", "mip-chance.log");
			env.start();
			
			// Create empty model
			GRBModel model = new GRBModel(env);			
			model.set(GRB.IntParam.LogToConsole, 0); // disable console logging
			
			// Create variables
		    GRBVar[] Q = new GRBVar[T];
		    GRBVar[] U = new GRBVar[T];
		    GRBVar[] z = new GRBVar[T];
		    GRBVar[] delta = new GRBVar[T];
		    for (int t = 0; t < T; t++) {
		    	Q[t] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");
		    	U[t] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "U");
		    	z[t] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "z");
		    	delta[t] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "delta");
		    }
		    		    
		    // expression variables
		    // staff flow
		    // x[t] = (x[t-1]+Q[t])(1-p)
		    GRBLinExpr[] x = new GRBLinExpr[T]; // end-of-period staff number
		    for (int t = 0; t < T; t++) {
		    	x[t] = new GRBLinExpr();
		    }
		    x[0].addTerm(1-turnoverRate[0], Q[0]);
		    for (int t = 1; t < T; t++) {
		    	x[t].multAdd(1-turnoverRate[t], x[t-1]);
		    	x[t].addTerm(1-turnoverRate[t], Q[t]);
		    }		    
		    
		    // objective function, set objective
		    GRBLinExpr obj = new GRBLinExpr();
		    for (int t = 0; t < T; t++) {
		    	obj.addTerm(fixCost, z[t]);
		    	obj.addTerm(unitVariCost, Q[t]);
		    	obj.addTerm(unitPenalty, U[t]);
		    	obj.multAdd(salary, x[t]);
		    }
		    model.setObjective(obj, GRB.MINIMIZE);
		    
		    // constraints
		    // Q_t <= z_t M
		    // U_t = (W-x_{t})^+:
		    // U_t <= w-x_{t}+(1-\delta)M
		    // U_t >= w-x_{t}-(1-\delta)M
		    // U_t <= \delta M
		    // w-x_{t} <= \delta M
		    int M = Integer.MAX_VALUE;
		    for (int t = 0; t < T; t++) {
		    	GRBLinExpr right1 = new GRBLinExpr();
		    	GRBLinExpr right2 = new GRBLinExpr();
		    	GRBLinExpr right3 = new GRBLinExpr();
		    	
		    	right1.addTerm(M, z[t]);
		    	model.addConstr(Q[t], GRB.LESS_EQUAL, right1, null); // Q_t <= z_t M
		    	
		    	GRBLinExpr deltaM = new GRBLinExpr();
		    	GRBLinExpr delta1M = new GRBLinExpr();
		    	GRBLinExpr deltaMinus1M = new GRBLinExpr();
		    	GRBLinExpr Wx = new GRBLinExpr();	
		    	
		    	deltaM.addTerm(M, delta[t]);
		    	delta1M.multAdd(-1, deltaM);
		    	delta1M.addConstant(M);
		    	deltaMinus1M.multAdd(1, deltaM);
		    	deltaMinus1M.addConstant(-M);
		    	Wx.addConstant(minStaffNum[t]);
		    	Wx.multAdd(-1, x[t]);
				
		    	right2.add(Wx);
		    	right2.add(delta1M);
		    	right3.add(Wx);
		    	right3.add(deltaMinus1M);
//		    	model.addConstr(U[t], GRB.GREATER_EQUAL, Wx, null); // may be only need this because of the objective
		    	model.addConstr(U[t], GRB.LESS_EQUAL, deltaM, null); // U_t <= delta M
		    	model.addConstr(U[t], GRB.LESS_EQUAL, right2, null); // U_t <= w-x_{t}+(1-\delta)M
		    	model.addConstr(U[t], GRB.GREATER_EQUAL, right3, null); // U_t >= w-x_{t-1}-Q-(1-\delta)M
		    	model.addConstr(Wx, GRB.LESS_EQUAL, deltaM, null); // w-x_{t-1}-Q <= \delta M	  
		    }
		    	
		    // Optimize model
			model.optimize();
			    
			// output results
			for (int t = 0; t < T; t++) {
				Z[t] = (int) z[t].get(GRB.DoubleAttr.X); // there is no specific get value method for gurobi in java     
		    }
			System.out.println(model.get(GRB.DoubleAttr.ObjVal));
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
		}
		
		return Z;
	}

}
