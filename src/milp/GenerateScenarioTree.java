package milp;

import java.util.Arrays;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBQuadExpr;
import gurobi.GRBVar;
import gurobi.GRB.IntParam;

/**
 * generate scenario tree for known distributions.
 * the objective is higher order function, guribo can not solve at present.
 * 
 * @author chen
 *
 */
public class GenerateScenarioTree {
	
	double[] QP(double mean, double std, int SampleNum) {
		try {
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv();	
			env.set(IntParam.OutputFlag, 0);
			env.start();
			
			// Create model
			GRBModel model = new GRBModel(env);
			model.set(GRB.IntParam.LogToConsole, 0);
			
			// Creat variables
			int N = SampleNum;
			GRBVar[] x = new GRBVar[N];
			GRBVar[] prob = new GRBVar[N];
			
			GRBQuadExpr meanExpec = new GRBQuadExpr();
			double[] ones = new double[N];
			Arrays.fill(ones, 1.0);
			meanExpec.addTerms(ones, x, prob);
			
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " +
			          e.getMessage());
		}		

		return null;
	}

	public static void main(String[] args) {
		double[] mean = {63, 27, 10, 24};
		double[] std = Arrays.stream(mean).map(i -> i*0.25).toArray();
		System.out.println(std[1]);

	}

}
