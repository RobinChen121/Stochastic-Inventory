package milp;


import org.apache.poi.hssf.record.LeftMarginRecord;
import org.apache.poi.poifs.crypt.DataSpaceMapUtils.IRMDSTransformInfo;

import cern.colt.Arrays;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBQConstr;
import gurobi.GRBQuadExpr;
import gurobi.GRBVar;

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
	
	public double pieceApprox(int[] z, int partionNum) {
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
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);		
			//env.set("logFile", "mip-chance.log");
			env.start();
			
			// Create empty model
			GRBModel model = new GRBModel(env);			
			//model.set(GRB.IntParam.LogToConsole, 0); // disable console logging
			
			// Create variables
		    GRBVar[] Q = new GRBVar[T];
		    GRBVar[] U = new GRBVar[T];
		    GRBVar[] xLb = new GRBVar[T]; 
		    //GRBVar[] delta = new GRBVar[T];
		    GRBVar[] x = new GRBVar[T];
		    for (int t = 0; t < T; t++) {
		    	Q[t] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");
		    	U[t] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "U");
		    	//delta[t] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "delta");
		    	x[t] = model.addVar(-Double.MAX_VALUE, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "x");
		    	xLb[t] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "y");
		    } 
		    
		   // objective function, set objective
		    GRBLinExpr obj = new GRBLinExpr();
		    for (int t = 0; t < T; t++) {
		    	if (z[t] == 1)
		    		obj.addConstant(fixCost);
		    	obj.addTerm(unitVariCost, Q[t]);
		    	obj.addTerm(unitPenalty, U[t]);
		    	obj.addTerm(salary, xLb[t]);
		    }
		    model.setObjective(obj, GRB.MINIMIZE);
		    
		    // z start periods
		    int[] zStart = new int[T];
		    for (int t = 0; t < T; t++) {
		    	if (z[t] == 1)
		    		zStart[t] = t;
		    	else {
		    		if (t == 0)
		    			zStart[t] = 0;
		    		else	
		    			zStart[t] = z[t-1];
				}
		    }
		    
		    // constraints
		    // Q_t <= z_t M
		    // U_t = (W-x_{t})^+:
		    // U_t >= W-x_{t}
		    // y = \sqrt((x[j-1]+Q[j])p(1-p))
		    int M = Integer.MAX_VALUE;
		    for (int t = 0; t < T; t++) {
		    	model.addConstr(Q[t], GRB.LESS_EQUAL, z[t]*M, null);// Q_t <= z_t M
		    	GRBLinExpr right1 = new GRBLinExpr();
		    	right1.addTerm(-1, x[t]);
		    	right1.addConstant(minStaffNum[t]);
		    		    	
		    	GRBLinExpr right2 = new GRBLinExpr();		    	
		    	int j = zStart[t];
		    	if (j == 0)
		    		right2.addTerm(turnoverRate[t]*(1-turnoverRate[t]), Q[j]);
		    	else {
		    		right2.addTerm(turnoverRate[t]*(1-turnoverRate[t]), Q[j]);
		    		right2.addTerm(turnoverRate[t]*(1-turnoverRate[t]), x[j-1]);
				}
		    	
		    	// piecewise constraints
		    	double sump = 0;
		    	double sumpmean = 0;
		    	for (int i = 0; i < partionNum; i++) {
		    		sump += prob[i];
		    		sumpmean += prob[i] * means[i];
		    		
		    		GRBLinExpr right1_1 = new GRBLinExpr();
		    	    right1_1.addTerm(sump, x[t]);
		    	    right1_1.addConstant(-4.8*sumpmean);
		    	    model.addConstr(xLb[t], GRB.GREATER_EQUAL, right1_1, null);		
		    	    
//		    	    GRBQuadExpr left1 = new GRBQuadExpr();
//			    	GRBQuadExpr left2 = new GRBQuadExpr();
//			    	GRBQuadExpr left3 = new GRBQuadExpr();
//			    	left1.addTerm(1, xLb[t], xLb[t]);
//			    	left2.addTerm(sump*sump, x[t], x[t]);
//			    	left3.addTerm(-2*sump, x[t], xLb[t]);
//			    	GRBQuadExpr left = new GRBQuadExpr();
//			    	left.add(left1);
//			    	left.add(left2);
//			    	left.add(left3);
//			    	
//			    	GRBLinExpr right = new GRBLinExpr();
//			    	right.multAdd(sumpmean*sumpmean, right2);
//			    	model.addQConstr(left, GRB.LESS_EQUAL, right, null);    
		    	}
		    	
		    	
		    	
		    	
		    	
		    	GRBLinExpr right3 = new GRBLinExpr();
		    	// U_t >= W-x_{t}
		    	right3.addTerm(-1, x[t]);
		    	right3.addConstant(minStaffNum[t]);
		    	model.addConstr(U[t], GRB.GREATER_EQUAL, right3, null); // U_t >= W-x_{t}
		    	
		    	// x_t+(x_{t-1}+Q_t)p-x_{t-1} >= 0
		    	// x_t+(x_{t-1}+Q_t)p-x_{t-1} <= z_tM
//		    	GRBLinExpr left4 = new GRBLinExpr();
//		    	left4.addTerm(turnoverRate[t], Q[t]);
//		    	if (t > 0) {
//					left4.addTerm(-(1-turnoverRate[t]), x[t-1]);
//				}
//		    	left4.addTerm(1, x[t]);
//		    	model.addConstr(left4, GRB.GREATER_EQUAL, 0, null);
//		    	model.addConstr(left4, GRB.LESS_EQUAL, z[t]*M, null); 
		    	
		    	// x_t = (x_{t-1}+Q_t)(1-p)
		    	GRBLinExpr right4 = new GRBLinExpr();
		    	right4.addTerm(1-turnoverRate[t], Q[t]);
		    	if (t > 0)
		    		right4.addTerm(1-turnoverRate[t], x[t-1]);
		    	//model.addConstr(xLb[t], GRB.GREATER_EQUAL, right4, null);
		    	model.addConstr(Q[t], GRB.EQUAL, 113, null);
		    	
		    }
		    
		    // Optimize model
			model.optimize();
			    
			// output results
			System.out.println(model.get(GRB.DoubleAttr.ObjVal));
			double[] Qx = new double[T];
			double[] xV = new double[T];
			double[] xLbV = new double[T];
			double[] Ux = new double[T];
			for (int t = 0; t < T; t++) {
				Qx[t] = Q[t].get(GRB.DoubleAttr.X);
				xV[t] = x[t].get(GRB.DoubleAttr.X);
				Ux[t] = U[t].get(GRB.DoubleAttr.X);
				xLbV[t] = xLb[t].get(GRB.DoubleAttr.X);
			}
			
			System.out.println("Q is " + Arrays.toString(Qx));
			System.out.println("x is " + Arrays.toString(xV));
			System.out.println("U is " + Arrays.toString(Ux));
			System.out.println("xLb is " + Arrays.toString(xLbV));
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
		}
		return 0;
	}

}
