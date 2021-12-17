package milp;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import sdp.cash.CashState;
import sdp.sampling.CartesianProduct;
import sdp.sampling.Sampling;
import umontreal.ssj.probdist.Distribution;

/**
 * @author chen
 * @email: okchen321@163.com
 * @date: 2021 Jun 24, 13:48:27  
 * @desp: calling gurobi in java to solve the chance-constrained problems/
 * 
 * the result is not good as expected, because the probability of non negative cash changes dramatically from
 * small to large.
 *
 */
public class PositiveCashChance {
	Distribution[] distributions;
	int[] sampleNums;
	double iniCash;
	double iniI;
	double[] price;
	double variCostUnit;
	double salvageValueUnit;
	double[] overheadCost;
	double serviceRate;
	double holdCostUnit;
	double[][] scenarios;
	
	double M1 = 10000;
	double M2 = 10000;
	
	double minB = 0;
	
	public PositiveCashChance(Distribution[] distributions, int[] sampleNums, double iniCash, double iniI, double[] price, double variCostUnit, 
			double salvageValueUnit, double holdCostUnit, double[] overheadCost, double serviceRate, double[][] scenarios) {
		this.distributions = distributions;
		this.sampleNums = sampleNums;
		this.iniCash = iniCash;
		this.iniI = iniI;
		this.price = price;
		this.variCostUnit = variCostUnit;
		this.salvageValueUnit = salvageValueUnit;
		this.overheadCost = overheadCost;
		this.serviceRate = serviceRate;
		this.scenarios = scenarios;
		this.holdCostUnit = holdCostUnit;
	}
	
	
	/**
	 * samples in each period form a scenario tree
	 * @return
	 */
	public double[] solve() {
		double[] result = new double[3];
		try {
			int T = distributions.length;
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			
			List<List<Integer>> Indexes = new ArrayList<>();
			for (int t = 0; t < T; t++)
				Indexes.add(IntStream.iterate(0, i -> i + 1).limit(sampleNums[t]).boxed().collect(Collectors.toList()));
			CartesianProduct<Integer> CP = new CartesianProduct<>();
			List<List<Integer>> scenarioIndexes = CP.product(Indexes);				
			
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);
			
			env.set("logFile", "mip-chance.log");
			env.start();

			// Create empty model
			GRBModel model = new GRBModel(env);
			
			model.set(GRB.IntParam.LogToConsole, 0);
			
			// Create variables
		    GRBVar[][] Q = new GRBVar[T][sampleNumTotal];		
		    GRBVar[][] delta = new GRBVar[T][sampleNumTotal];	// auxiliary variable
		    GRBVar[][] I = new GRBVar[T][sampleNumTotal];	
		    for (int t = 0; t < T; t++) {
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		Q[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");
		    		delta[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "delta");
		    		I[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "I");		    		
		    	}
		    }		    
		    
		    // expression variables
		    GRBLinExpr[][] cash = new GRBLinExpr[T][sampleNumTotal];
		    GRBLinExpr[][] minCash = new GRBLinExpr[T][sampleNumTotal];
		    GRBLinExpr[][] revenue = new GRBLinExpr[T][sampleNumTotal];
		    GRBLinExpr expectFinalCash = new GRBLinExpr();
		    
		    // cash flow
		    for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		revenue[t][s] = new GRBLinExpr();
		    		cash[t][s] = new GRBLinExpr();
		    		minCash[t][s] = new GRBLinExpr();
		    		if (t == 0 && T > 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}
		    		else if (t == 0 && T == 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else if (t == T - 1 && T > 1) {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);	
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}	    			
		    	}
		    
		    // objective function, set objective
		    expectFinalCash.clear();
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	expectFinalCash.multAdd(1.0/sampleNumTotal, cash[T - 1][s]);
		    }
			model.setObjective(expectFinalCash, GRB.MAXIMIZE);

			// Add constraint
			// inventory flow
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		int sIndex = scenarioIndexes.get(s).get(t);
		    		double demand = scenarios[t][sIndex]; 
		    		if (t == 0) {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
		    			rightExpr1.addConstant(iniI); rightExpr1.addTerm(1, Q[t][s]); 
		    			rightExpr1.addConstant(-demand); 
		    			GRBLinExpr rightExpr2 = new GRBLinExpr();
		    			rightExpr2.addConstant(M1); rightExpr2.addTerm(-M1, delta[t][s]);
		    			GRBLinExpr rightExpr = new GRBLinExpr();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr2);
		    			model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint1");
		    			
		    			GRBLinExpr rightExpr3 = new GRBLinExpr();
		    			rightExpr3.multAdd(-1, rightExpr2);
		    			rightExpr.clear();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr3);
		    			model.addConstr(I[t][s], GRB.GREATER_EQUAL, rightExpr, "IConstraint2");
		    			
		    			GRBLinExpr rightExpr4 = new GRBLinExpr();
		    			rightExpr4.addTerm(M1, delta[t][s]);
		    			model.addConstr(rightExpr1, GRB.LESS_EQUAL, rightExpr4, "IConstraint3");
		    		}
		    		else {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
		    			rightExpr1.addTerm(1, I[t-1][s]); rightExpr1.addTerm(1, Q[t][s]); 
		    			rightExpr1.addConstant(-demand); 
		    			GRBLinExpr rightExpr2 = new GRBLinExpr();
		    			rightExpr2.addConstant(M1); rightExpr2.addTerm(-M1, delta[t][s]);
		    			GRBLinExpr rightExpr = new GRBLinExpr();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr2);
		    			model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint1");
		    			
		    			GRBLinExpr rightExpr3 = new GRBLinExpr();
		    			rightExpr3.multAdd(-1, rightExpr2);
		    			rightExpr.clear();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr3);
		    			model.addConstr(I[t][s], GRB.GREATER_EQUAL, rightExpr, "IConstraint2");
		    			
		    			GRBLinExpr rightExpr4 = new GRBLinExpr();
		    			rightExpr4.addTerm(M1, delta[t][s]);
		    			model.addConstr(rightExpr1, GRB.LESS_EQUAL, rightExpr4, "IConstraint3");
					}
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(M1, delta[t][s]);
		    		model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint4");
		    	}
			
			// chance constraint
		    GRBVar[] alpha = new GRBVar[sampleNumTotal];	// whether cash balance is negative in this scenario
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	alpha[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "alpha");
		    }
			GRBLinExpr sumAlpha = new GRBLinExpr();
			sumAlpha.clear();
			for (int s = 0; s < sampleNumTotal; s++) {
				sumAlpha.addTerm(1, alpha[s]);
			}
			model.addConstr(sumAlpha, GRB.LESS_EQUAL, negativeScenarioNumRequire, "ChanceConstraint1");
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
//		    		if(t == 0) {
//		    			minCash[t][s].addConstant(iniCash); 
//		    			minCash[t][s].addTerm(-variCostUnit, Q[t][s]);
//		    			// minCash[t][s].addConstant(-overheadCost[t]);
//		    		}
//		    		else {
//		    			minCash[t][s].multAdd(1, cash[t-1][s]); 
//		    			minCash[t][s].addTerm(-variCostUnit, Q[t][s]);
////		    			// minCash[t][s].addConstant(-overheadCost[t]);
//					}
		    		minCash[t][s].multAdd(1, cash[t][s]); 
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(-M2, alpha[s]);
		    		rightExpr.addConstant(minB);
		    		model.addConstr(minCash[t][s], GRB.GREATER_EQUAL, rightExpr, "ChanceConstraint2");
		    	}
			
			// first-stage decision, here and now decision
			for (int s = 0; s < sampleNumTotal - 1; s++) {
				model.addConstr(Q[0][s], GRB.EQUAL, Q[0][s+1], "fistStageConstraint");
			}
			
			// Optimize model
		    model.optimize();
		    
		    int negScenarioNum = 0;
		    double[] alphaValue = new double[sampleNumTotal];
		    
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter("detail-results.txt"));
				out.write("ordering quantity and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", Q[t][s].get(GRB.DoubleAttr.X)));
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("cash balance in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", cash[t][s].getValue()));
					alphaValue[s] = alpha[s].get(GRB.DoubleAttr.X);
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("min cash balance in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					int recordBefore = 0;
					for (int t = 0; t < T; t++) {
						out.write(String.format("%.1f  ", minCash[t][s].getValue()));
						if (recordBefore == 0 && minCash[t][s].getValue() < -0.1 ){ // a minus value to count correctly
							negScenarioNum++;
							recordBefore = 1;
						}
					}
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("Inventory balance in each period and each scenario: \n");
			    for (int s = 0; s < sampleNumTotal; s++) {
			    	out.write("scenario " + s + ": \n");
			    	for (int t = 0; t < T; t++)
			    		out.write(String.format("%.1f  ", I[t][s].get(GRB.DoubleAttr.X)));
			    	out.newLine();
			    }
				out.close();
			} catch (IOException e) {				
			}	    
		    
//			System.out.println();
//			System.out.println("negative scenario number in the solution is : " + negScenarioNum);
//			System.out.println("maximum negative scenario number required is: " + negativeScenarioNumRequire);
//			System.out.println("total scenario number is : " + sampleNumTotal);
			
			result[0] = Q[0][0].get(GRB.DoubleAttr.X);
			result[1] = model.get(GRB.DoubleAttr.ObjVal);
			result[2] = negScenarioNum;
			
			// Dispose of model and environment
			model.dispose();
			env.dispose();
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			result[0] = 0;
			result[1] = 0;
			result[2] = 0;
		}
		return result;
	}
	
	/**
	 * solve the chance SAA by sorting the scenarios;
	 * objective is to maximize expected final cash
	 * @return
	 */
	public double[] solveSort() {
		double[] result = new double[3];
		try {
			int T = distributions.length;
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			
			List<List<Double>> listScenarios = new ArrayList<>();
			for (int t = 0; t < T; t++)
				listScenarios.add(Arrays.stream(scenarios[t]).boxed().collect(Collectors.toList()));
			CartesianProduct<Double> CP = new CartesianProduct<>();
			List<List<Double>> listScenariosAll = CP.product(listScenarios);
			Comparator<List<Double>> comparator = (o1, o2) -> {
				double sum1 = 0;
				double sum2 = 0;
				for (int t = 0; t < o1.size(); t++) {
					sum1 += o1.get(t)*price[t];
					sum2 += o2.get(t)*price[t];
				}
				if (sum1 < sum2)
					return 1;
				else if(sum1 > sum2)
					return -1;
				else {
					return 0;
				}				
			};
			Collections.sort(listScenariosAll, comparator);
			
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);
			env.set("logFile", "mip-chance.log");
			env.start();

			// Create empty model
			GRBModel model = new GRBModel(env);
			model.set(GRB.IntParam.LogToConsole, 0);

			// Create variables
			GRBVar[][] Q = new GRBVar[T][sampleNumTotal];
			GRBVar[][] delta = new GRBVar[T][sampleNumTotal]; // auxiliary variable
			GRBVar[][] I = new GRBVar[T][sampleNumTotal];
			for (int t = 0; t < T; t++) {
				for (int s = 0; s < sampleNumTotal; s++) {
					Q[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");
					delta[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "delta");
					I[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "I");
				}
			}

			// expression variables
			GRBLinExpr[][] cash = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr[][] minCash = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr[][] revenue = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr expectFinalCash = new GRBLinExpr();		
			
			// cash flow
		    for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		revenue[t][s] = new GRBLinExpr();
		    		cash[t][s] = new GRBLinExpr();
		    		minCash[t][s] = new GRBLinExpr();
		    		if (t == 0 && T > 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}
		    		else if (t == 0 && T == 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else if (t == T - 1 && T > 1) {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);	
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}	    			
		    	}
		    
		    // objective function, set objective
		    expectFinalCash.clear();
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	expectFinalCash.multAdd(1.0/sampleNumTotal, cash[T - 1][s]);
		    }
			model.setObjective(expectFinalCash, GRB.MAXIMIZE);

			// Add constraint
			// inventory flow
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		double demand = listScenariosAll.get(s).get(t);
		    		if (t == 0) {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
		    			rightExpr1.addConstant(iniI); rightExpr1.addTerm(1, Q[t][s]); 
		    			rightExpr1.addConstant(-demand); 
		    			GRBLinExpr rightExpr2 = new GRBLinExpr();
		    			rightExpr2.addConstant(M1); rightExpr2.addTerm(-M1, delta[t][s]);
		    			GRBLinExpr rightExpr = new GRBLinExpr();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr2);
		    			model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint1");
		    			
		    			GRBLinExpr rightExpr3 = new GRBLinExpr();
		    			rightExpr3.multAdd(-1, rightExpr2);
		    			rightExpr.clear();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr3);
		    			model.addConstr(I[t][s], GRB.GREATER_EQUAL, rightExpr, "IConstraint2");
		    			
		    			GRBLinExpr rightExpr4 = new GRBLinExpr();
		    			rightExpr4.addTerm(M1, delta[t][s]);
		    			model.addConstr(rightExpr1, GRB.LESS_EQUAL, rightExpr4, "IConstraint3");
		    		}
		    		else {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
		    			rightExpr1.addTerm(1, I[t-1][s]); rightExpr1.addTerm(1, Q[t][s]); 
		    			rightExpr1.addConstant(-demand); 
		    			GRBLinExpr rightExpr2 = new GRBLinExpr();
		    			rightExpr2.addConstant(M1); rightExpr2.addTerm(-M1, delta[t][s]);
		    			GRBLinExpr rightExpr = new GRBLinExpr();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr2);
		    			model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint1");
		    			
		    			GRBLinExpr rightExpr3 = new GRBLinExpr();
		    			rightExpr3.multAdd(-1, rightExpr2);
		    			rightExpr.clear();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr3);
		    			model.addConstr(I[t][s], GRB.GREATER_EQUAL, rightExpr, "IConstraint2");
		    			
		    			GRBLinExpr rightExpr4 = new GRBLinExpr();
		    			rightExpr4.addTerm(M1, delta[t][s]);
		    			model.addConstr(rightExpr1, GRB.LESS_EQUAL, rightExpr4, "IConstraint3");
					}
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(M1, delta[t][s]);
		    		model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint4");
		    	}
			
			// chance constraint
//		    GRBVar[] alpha = new GRBVar[sampleNumTotal];	// whether cash balance is negative in this scenario
//		    for (int s = 0; s < sampleNumTotal; s++) {
//		    	alpha[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "alpha");
//		    }
//			GRBLinExpr sumAlpha = new GRBLinExpr();
//			sumAlpha.clear();
//			for (int s = 0; s < sampleNumTotal; s++) {
//				sumAlpha.addTerm(1, alpha[s]);
//			}
//			model.addConstr(sumAlpha, GRB.LESS_EQUAL, negativeScenarioNumRequire, "ChanceConstraint1");
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
//		    		if(t == 0) {
//	    				minCash[t][s].addConstant(iniCash); 
//	    				minCash[t][s].addTerm(-variCostUnit, Q[t][s]);
//	    				// minCash[t][s].addConstant(-overheadCost[t]);
//	    			}
//	    			else {
//	    				minCash[t][s].multAdd(1, cash[t-1][s]); 
//	    				minCash[t][s].addTerm(-variCostUnit, Q[t][s]);
////	    			// minCash[t][s].addConstant(-overheadCost[t]);
//					}
		    		minCash[t][s].multAdd(1, cash[t][s]); 
		    		
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		if (s < negativeScenarioNumRequire) { // no =
		    			rightExpr.addConstant(-M2);
		    			rightExpr.addConstant(minB);
		    		}
		    			
		    		else {
						rightExpr.clear();
					}
		    		model.addConstr(minCash[t][s], GRB.GREATER_EQUAL, rightExpr, "ChanceConstraint2");
		    	}
			
			// first-stage decision, here and now decision
			for (int s = 0; s < sampleNumTotal - 1; s++) {
				model.addConstr(Q[0][s], GRB.EQUAL, Q[0][s+1], "fistStageConstraint");
			}
			
			// Optimize model
		    model.optimize();
		    
		    int negScenarioNum = 0;
		    double[] alphaValue = new double[sampleNumTotal];
		    
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter("detail-results.txt"));
				out.write("ordering quantity and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", Q[t][s].get(GRB.DoubleAttr.X)));
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("cash balance in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", cash[t][s].getValue()));
//					alphaValue[s] = alpha[s].get(GRB.DoubleAttr.X);
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("min cash balance in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					int recordBefore = 0;
					for (int t = 0; t < T; t++) {
						out.write(String.format("%.1f  ", minCash[t][s].getValue()));
						if (recordBefore == 0 && minCash[t][s].getValue() < -0.1 ){ // a minus value to count correctly
							negScenarioNum++;
							recordBefore = 1;
						}
					}
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("Inventory balance in each period and each scenario: \n");
			    for (int s = 0; s < sampleNumTotal; s++) {
			    	out.write("scenario " + s + ": \n");
			    	for (int t = 0; t < T; t++)
			    		out.write(String.format("%.1f  ", I[t][s].get(GRB.DoubleAttr.X)));
			    	out.newLine();
			    }
				out.close();
			} catch (IOException e) {				
			}	    
		    
//			System.out.println();
//			System.out.println("negative scenario number in the solution is : " + negScenarioNum);
//			System.out.println("maximum negative scenario number required is: " + negativeScenarioNumRequire);
//			System.out.println("total scenario number is : " + sampleNumTotal);
			
			result[0] = Q[0][0].get(GRB.DoubleAttr.X);
			result[1] = model.get(GRB.DoubleAttr.ObjVal);
			result[2] = negScenarioNum;
			
			// Dispose of model and environment
			model.dispose();
			env.dispose();
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			result[0] = 0;
			result[1] = 0;
			result[2] = 0;
		}
		return result;
	}
	
	
	/**
	 * 
	 * the samples are all independent, do not form a tree.
	 * samples in all periods are sampled once all;
	 * running slow sometimes
	 * @return
	 */
	public double[] solve2(int fixSampleNum) {
		double[] result = new double[3];
		try {
			int T = distributions.length;
			int sampleNumTotal = fixSampleNum;
			Sampling sampling = new Sampling();
			double[][] scenarios = sampling.generateLHSamples(distributions, sampleNumTotal); 
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));			
			
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);
			env.set("logFile", "mip-chance.log");
			env.start();

			// Create empty model
			GRBModel model = new GRBModel(env);

			// Create variables
		    GRBVar[][] Q = new GRBVar[T][sampleNumTotal];		
		    GRBVar[][] delta = new GRBVar[T][sampleNumTotal];	// auxiliary variable
		    GRBVar[][] I = new GRBVar[T][sampleNumTotal];	
		    for (int t = 0; t < T; t++) {
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		Q[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");
		    		delta[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "delta");
		    		I[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "I");		    		
		    	}
		    }		    
		    
		    // expression variables
		    GRBLinExpr[][] cash = new GRBLinExpr[T][sampleNumTotal];
		    GRBLinExpr[][] minCash = new GRBLinExpr[T][sampleNumTotal];
		    GRBLinExpr[][] revenue = new GRBLinExpr[T][sampleNumTotal];
		    GRBLinExpr expectFinalCash = new GRBLinExpr();
		    
		    // cash flow
		    for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		revenue[t][s] = new GRBLinExpr();
		    		cash[t][s] = new GRBLinExpr();
		    		minCash[t][s] = new GRBLinExpr();
		    		if (t == 0 && T > 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}
		    		else if (t == 0 && T == 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else if (t == T - 1 && T > 1) {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit, Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);	
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}	    			
		    	}
		    
		    // objective function, set objective
		    expectFinalCash.clear();
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	expectFinalCash.multAdd(1.0/sampleNumTotal, cash[T - 1][s]);
		    }
			model.setObjective(expectFinalCash, GRB.MAXIMIZE);

			// Add constraint
			// inventory flow
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		double demand = scenarios[s][t];
		    		if (t == 0) {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
		    			rightExpr1.addConstant(iniI); rightExpr1.addTerm(1, Q[t][s]); 
		    			rightExpr1.addConstant(-demand); 
		    			GRBLinExpr rightExpr2 = new GRBLinExpr();
		    			rightExpr2.addConstant(M1); rightExpr2.addTerm(-M1, delta[t][s]);
		    			GRBLinExpr rightExpr = new GRBLinExpr();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr2);
		    			model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint1");
		    			
		    			GRBLinExpr rightExpr3 = new GRBLinExpr();
		    			rightExpr3.multAdd(-1, rightExpr2);
		    			rightExpr.clear();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr3);
		    			model.addConstr(I[t][s], GRB.GREATER_EQUAL, rightExpr, "IConstraint2");
		    			
		    			GRBLinExpr rightExpr4 = new GRBLinExpr();
		    			rightExpr4.addTerm(M1, delta[t][s]);
		    			model.addConstr(rightExpr1, GRB.LESS_EQUAL, rightExpr4, "IConstraint3");
		    		}
		    		else {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
		    			rightExpr1.addTerm(1, I[t-1][s]); rightExpr1.addTerm(1, Q[t][s]); 
		    			rightExpr1.addConstant(-demand); 
		    			GRBLinExpr rightExpr2 = new GRBLinExpr();
		    			rightExpr2.addConstant(M1); rightExpr2.addTerm(-M1, delta[t][s]);
		    			GRBLinExpr rightExpr = new GRBLinExpr();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr2);
		    			model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint1");
		    			
		    			GRBLinExpr rightExpr3 = new GRBLinExpr();
		    			rightExpr3.multAdd(-1, rightExpr2);
		    			rightExpr.clear();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr3);
		    			model.addConstr(I[t][s], GRB.GREATER_EQUAL, rightExpr, "IConstraint2");
		    			
		    			GRBLinExpr rightExpr4 = new GRBLinExpr();
		    			rightExpr4.addTerm(M1, delta[t][s]);
		    			model.addConstr(rightExpr1, GRB.LESS_EQUAL, rightExpr4, "IConstraint3");
					}
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(M1, delta[t][s]);
		    		model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint4");
		    	}
			
			// chance constraint
		    GRBVar[] alpha = new GRBVar[sampleNumTotal];	// whether cash balance is negative in this scenario
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	alpha[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "alpha");
		    }
			GRBLinExpr sumAlpha = new GRBLinExpr();
			sumAlpha.clear();
			for (int s = 0; s < sampleNumTotal; s++) {
				sumAlpha.addTerm(1, alpha[s]);
			}
			model.addConstr(sumAlpha, GRB.LESS_EQUAL, negativeScenarioNumRequire, "ChanceConstraint1");
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		if(t == 0) {
		    			minCash[t][s].addConstant(iniCash); minCash[t][s].addTerm(-variCostUnit, Q[t][s]);
		    			// minCash[t][s].addConstant(-overheadCost[t]);
		    		}
		    		else {
		    			minCash[t][s].multAdd(1, cash[t-1][s]); minCash[t][s].addTerm(-variCostUnit, Q[t][s]);
		    			// minCash[t][s].addConstant(-overheadCost[t]);
					}
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(-M2, alpha[s]);
		    		model.addConstr(minCash[t][s], GRB.GREATER_EQUAL, rightExpr, "ChanceConstraint2");
		    	}
			
			// first-stage decision, here and now decision
			for (int s = 0; s < sampleNumTotal - 1; s++) {
				model.addConstr(Q[0][s], GRB.EQUAL, Q[0][s+1], "fistStageConstraint");
			}
			
			// Optimize model
		    model.optimize();
		    
		    int negScenarioNum = 0;
		    double[] alphaValue = new double[sampleNumTotal];
		    
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter("detail-results2.txt"));
				out.write("ordering quantity and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", Q[t][s].get(GRB.DoubleAttr.X)));
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("cash balance in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", cash[t][s].getValue()));
					alphaValue[s] = alpha[s].get(GRB.DoubleAttr.X);
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("min cash balance in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					int recordBefore = 0;
					for (int t = 0; t < T; t++) {
						out.write(String.format("%.1f  ", minCash[t][s].getValue()));
						if (recordBefore == 0 && minCash[t][s].getValue() < -0.1 ){ // a minus value to count correctly
							negScenarioNum++;
							recordBefore = 1;
						}
					}
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("Inventory balance in each period and each scenario: \n");
			    for (int s = 0; s < sampleNumTotal; s++) {
			    	out.write("scenario " + s + ": \n");
			    	for (int t = 0; t < T; t++)
			    		out.write(String.format("%.1f  ", I[t][s].get(GRB.DoubleAttr.X)));
			    	out.newLine();
			    }
				out.close();
			} catch (IOException e) {				
			}	    
		    
			System.out.println();
			System.out.println("negative scenario number in the solution is : " + negScenarioNum);
			System.out.println("maximum negative scenario number required is: " + negativeScenarioNumRequire);
			System.out.println("total scenario number is : " + sampleNumTotal);
			
			result[0] = Q[0][0].get(GRB.DoubleAttr.X);
			result[1] = model.get(GRB.DoubleAttr.ObjVal);
			result[2] = negScenarioNum;
			
			// Dispose of model and environment
			model.dispose();
			env.dispose();
		} catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			result[0] = 0;
			result[1] = 0;
			result[2] = 0;
		}
		return result;
	}

}
