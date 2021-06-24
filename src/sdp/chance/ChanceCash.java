package sdp.chance;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.IntPredicate;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import milp.GurobiChance;
import sdp.sampling.CartesianProduct;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;


/**
 * @author chen
 * @email: okchen321@163.com
 * @date: 2021 Jun 18, 17:41:22  
 * @desp: try the chance programming approach for the stochastic inventory problem.
 *
 */
public class ChanceCash {
	
	static int[][] scenarioIndexes(int[] sampleNums, int sampleNumTotal){
		int T = sampleNums.length;
		int[][] arr = new int[sampleNumTotal][T];
		for (int i = 0; i < sampleNumTotal; i++) {
			for (int t = 0; t < T; t++) {
				arr[i][t] = 0;
			}
		}
		return null;
	}
	
	public static void main(String[] args) {
		double iniCash = 100;
		double iniI = 0;
		double[] price = {6, 3, 5};
		double variCostUnit = 2;
		double salvageValueUnit = 0.5;
		double trunQuantile = 0.9999;
		double serviceRate = 0.9; // maximum negative possibility rate is 1 - serviceRate

		double[] meanDemand = {40, 30, 30};
		int[] sampleNums = {10, 10, 10}; // sample number in each period
		double maxOrderQuantity = 200; // maximum ordering quantity when having enough cash
		int T = sampleNums.length;
		double[] overheadCost = new double[T];
		Arrays.fill(overheadCost, 100); // overhead costs
		

		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				// .mapToObj(i -> new NormalDist(meanDemand[i], Math.sqrt(meanDemand[i]))) //
				// can be changed to other distributions
				.mapToObj(i -> new PoissonDist(meanDemand[i])).toArray(Distribution[]::new);

		
		GurobiChance model = new GurobiChance(distributions, sampleNums, iniCash, iniI, price, variCostUnit, 
							salvageValueUnit, overheadCost, serviceRate);
		long currTime = System.currentTimeMillis();
		double[] result = model.solve();
		System.out.println();
	    System.out.println("**********************************************");
	    System.out.println("result of SAA: ");
		double time = (System.currentTimeMillis() - currTime) / 1000.00;
		System.out.println("running time is " + time + "s");	
	    System.out.printf("first stage decison Q is: %.2f\n", result[0]);
	    System.out.printf("Objective value is: %.2f\n", result[1]);
	}
}
