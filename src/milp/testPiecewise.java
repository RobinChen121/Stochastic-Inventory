package milp;

import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

/**
* @author Zhen Chen
* @date: 2018-11-16, pm 4:53:57  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  test the usage of piecewise function of cplex in java
*/

public class testPiecewise {

	public static void main(String[] args) {
		try {
			IloCplex cplex = new IloCplex();
			
			IloNumVar x = cplex.numVar(-Double.MAX_VALUE, Double.MAX_VALUE);
			
			double[] points = {100, 200};
			double[] slopes = {1, 2, -3};
			IloNumExpr fx = cplex.piecewiseLinear(x, points, slopes, 0, 300);
						
			cplex.addMaximize(fx);
			
			cplex.addEq(x, 200);
			
			if (cplex.solve()) {
				System.out.println(cplex.getValue(x));
				System.out.println("Solution value = " + cplex.getObjValue());
			}
			
		} catch (Exception e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		

	}

}


