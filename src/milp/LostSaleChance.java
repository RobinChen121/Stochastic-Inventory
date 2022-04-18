package milp;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.IntPredicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;

import cern.jet.random.Beta;
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
 * @desp: calling gurobi in java to solve the chance-constrained problems
 * 
 * the lost sale rate is a joint chance constraint. 
 *
 */
public class LostSaleChance {
	Distribution[] distributions;
	int[] sampleNums;
	double iniCash;
	double iniI;
	double[] price;
	double[] variCostUnit;
	double salvageValueUnit;
	double[] overheadCost;
	double serviceRate;
	double holdCostUnit;
	double[][] scenarios;
	
	double M1; // revise
	double M2;
	double M3;
	
	
	public LostSaleChance(Distribution[] distributions, int[] sampleNums, double iniCash, double iniI, double[] price, double[] variCostUnit, 
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
	 * initial SAA model for the joint chance problem.
	 * samples in each period form a scenario tree.
	 * @return survival probability
	 */
	public double[] solveMaxSurvival() {
		double[] result = new double[3];
		try {
			int T = distributions.length;
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			
			List<List<Integer>> Indexes = new ArrayList<>();
			for (int t = 0; t < T; t++)
				Indexes.add(IntStream.iterate(0, i -> i + 1).limit(sampleNums[t]).boxed().collect(Collectors.toList()));
			CartesianProduct<Integer> CP = new CartesianProduct<>();
			List<List<Integer>> scenarioIndexes = CP.product(Indexes); // too slow when T >= 10		
			
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);
			
			//env.set("logFile", "mip-chance.log");
			env.start();

			// Create empty model
			GRBModel model = new GRBModel(env);
			
			model.set(GRB.IntParam.LogToConsole, 0);
			
			// Create variables
		    GRBVar[][] Q = new GRBVar[T][sampleNumTotal];		
		    GRBVar[][] delta = new GRBVar[T][sampleNumTotal];	// auxiliary variable
		    GRBVar[][] I = new GRBVar[T][sampleNumTotal];	
		    GRBVar[][] alpha = new GRBVar[T][sampleNumTotal];	// auxiliary variable
		    GRBVar[] z = new GRBVar[sampleNumTotal];	// auxiliary variable
		    GRBVar[] beta = new GRBVar[sampleNumTotal];	// whether lost sale happens
		    
		    
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	for (int t = 0; t < T; t++) {
		    		Q[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");
		    		delta[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "delta");
		    		I[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "I");	
		    		alpha[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "alpha");
		    	}
		    	z[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "z");
		    	beta[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "beta");
		    }		    
		    
		    // expression variables
		    GRBLinExpr[][] cash = new GRBLinExpr[T][sampleNumTotal];
		    GRBLinExpr[][] revenue = new GRBLinExpr[T][sampleNumTotal];
		    GRBLinExpr expectFinalPositiveCashNum = new GRBLinExpr();
		    
		    // cash flow
			// choose maximum sum of demand in T periods as M1
			M1 = 0;
			for (int s = 0; s < sampleNumTotal; s++) {
				double thissSumD = 0;
				for (int t = 0; t < T; t++) {
					int sIndex = scenarioIndexes.get(s).get(t);
					thissSumD += scenarios[t][sIndex];
				}
				if (thissSumD > M1)
					M1 = thissSumD;
			}
			M2 = iniCash + price[0] * M1; // prices are same
			M3 = iniI * holdCostUnit * T + variCostUnit[0] * M1 + Arrays.stream(overheadCost).sum() - iniCash; // variCostUnit[0] * M1
		    for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		revenue[t][s] = new GRBLinExpr();
		    		cash[t][s] = new GRBLinExpr();
		    		if (t == 0 && T > 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); 
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    			cash[t][s].addConstant(-overheadCost[t]);
		    		}
		    		else if (t == 0 && T == 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    			 cash[t][s].addConstant(-overheadCost[t]);
					}
		    		else if (t == T - 1 && T > 1) {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); 
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    			cash[t][s].addConstant(-overheadCost[t]);
					}
		    		else {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); 	
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    			cash[t][s].addConstant(-overheadCost[t]);
		    		}	    			
		    	}
		    
		    // objective function, set objective
		    expectFinalPositiveCashNum.clear();
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	expectFinalPositiveCashNum.addTerm(1.0, z[s]); // make objective greater 
		    }
			model.setObjective(expectFinalPositiveCashNum, GRB.MAXIMIZE);
			
			

			//M1 = M1 + iniI;
			
			// Add constraint
			// inventory flow
			for (int t = 0; t < T; t++) {	
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		int sIndex = scenarioIndexes.get(s).get(t);
		    		double demand = scenarios[t][sIndex]; 
		    		if (t == 0) {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
		    			rightExpr1.addConstant(iniI); rightExpr1.addTerm(1, Q[t][s]); 
		    			rightExpr1.addConstant(-demand); 
		    			GRBLinExpr rightExpr2 = new GRBLinExpr();
		    			rightExpr2.addTerm(M1, delta[t][s]);
		    			GRBLinExpr rightExpr = new GRBLinExpr();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr2);
		    			model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint1");
		    			
		    			GRBLinExpr rightExpr3 = new GRBLinExpr();
		    			rightExpr3.multAdd(-1, rightExpr2);
		    			rightExpr.clear();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr3);
		    			model.addConstr(I[t][s], GRB.GREATER_EQUAL, rightExpr, "IConstraint2");
		    			
		    			GRBLinExpr rightExpr4 = new GRBLinExpr();
		    			rightExpr4.addConstant(M1); rightExpr4.addTerm(-M1, delta[t][s]);
		    			model.addConstr(rightExpr1, GRB.LESS_EQUAL, rightExpr4, "IConstraint3");
		    		}
		    		else {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
		    			rightExpr1.addTerm(1, I[t-1][s]); rightExpr1.addTerm(1, Q[t][s]); 
		    			rightExpr1.addConstant(-demand); 
		    			GRBLinExpr rightExpr2 = new GRBLinExpr();
		    			rightExpr2.addTerm(M1, delta[t][s]);
		    			GRBLinExpr rightExpr = new GRBLinExpr();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr2);
		    			model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint1");
		    			
		    			GRBLinExpr rightExpr3 = new GRBLinExpr();
		    			rightExpr3.multAdd(-1, rightExpr2);
		    			rightExpr.clear();
		    			rightExpr.add(rightExpr1); rightExpr.add(rightExpr3);
		    			model.addConstr(I[t][s], GRB.GREATER_EQUAL, rightExpr, "IConstraint2");
		    			
		    			GRBLinExpr rightExpr4 = new GRBLinExpr();
		    			rightExpr4.addConstant(M1); rightExpr4.addTerm(-M1, delta[t][s]);
		    			model.addConstr(rightExpr1, GRB.LESS_EQUAL, rightExpr4, "IConstraint3");
					}
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addConstant(M1); rightExpr.addTerm(-M1, delta[t][s]);
		    		model.addConstr(I[t][s], GRB.LESS_EQUAL, rightExpr, "IConstraint4");
		    	}
			}
			
			// chance constraint
			for (int s = 0; s < sampleNumTotal; s++) {
				for (int t = 0; t < T; t++) {
					model.addConstr(delta[t][s], GRB.LESS_EQUAL, beta[s], "ChanceConstraint1");
				}
			}
			
			GRBLinExpr sumBeta = new GRBLinExpr();	
			for (int s = 0; s < sampleNumTotal; s++) {
				sumBeta.addTerm(1, beta[s]);
			}
			model.addConstr(sumBeta, GRB.LESS_EQUAL, negativeScenarioNumRequire, "ChanceConstraint2");
			

			
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++){
//		    		model.addGenConstrIndicator(alpha[t][s], 1, cash[t][s], GRB.GREATER_EQUAL, 0, null);
//		    		model.addGenConstrIndicator(alpha[t][s], 0, cash[t][s], GRB.LESS_EQUAL, 0, null);
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(M2, alpha[t][s]); rightExpr.addConstant(-2); // add a constant to avoid survival probs too big
		    		model.addConstr(cash[t][s], GRB.LESS_EQUAL, rightExpr, "CashConstraint1");
		    		rightExpr.clear();
		    		rightExpr.addConstant(-M3); rightExpr.addTerm(M3, alpha[t][s]);
		    		model.addConstr(cash[t][s], GRB.GREATER_EQUAL, rightExpr, "CashConstraint4");
		    	}
			
			for (int s = 0; s < sampleNumTotal; s++) {
				GRBLinExpr leftExpr = new GRBLinExpr();
    			leftExpr.addTerm(-1, z[s]); leftExpr.addConstant(1);
    			GRBLinExpr rightExpr = new GRBLinExpr();
    			//rightExpr.clear();
		    	for (int t = 0; t < T; t++){		    
		    		model.addConstr(z[s], GRB.LESS_EQUAL, alpha[t][s], "CashConstraint3");
		    		rightExpr.addTerm(-1, alpha[t][s]);  rightExpr.addConstant(1);
		    	}
		    	model.addConstr(leftExpr, GRB.LESS_EQUAL, rightExpr, "CashConstraint2");
			}
			
			// first-stage decision, here and now decision
			for (int s = 0; s < sampleNumTotal - 1; s++) {
				model.addConstr(Q[0][s], GRB.EQUAL, Q[0][s+1], "fistStageConstraint");
			}
			
			// Optimize model
			//model.set(GRB.DoubleParam.TimeLimit, 3600); // set time limit for gurobi
		    model.optimize();
		    
		    int negScenarioNum = 0;
		    double[] zValue = new double[sampleNumTotal];
		    double[] betaValue = new double[sampleNumTotal];
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
					int thissNegative = 0; 
					for (int t = 0; t < T; t++) {
						out.write(String.format("%.1f  ", cash[t][s].getValue()));
						if (cash[t][s].getValue() < -1)
							thissNegative = 1;
					}
					negScenarioNum += thissNegative;
					zValue[s] = z[s].get(GRB.DoubleAttr.X);
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("alpha in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", alpha[t][s].get(GRB.DoubleAttr.X)));
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("z in each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write(String.format("%.1f  ", z[s].get(GRB.DoubleAttr.X)));	
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
			    	betaValue[s] = beta[s].get(GRB.DoubleAttr.X);
			    }
				out.close();
			} catch (IOException e) {				
			}	    
			
			result[0] = Q[0][0].get(GRB.DoubleAttr.X);
			result[1] = model.get(GRB.DoubleAttr.ObjVal);
			result[2] = Arrays.stream(betaValue).sum();
			//System.out.println("negative cash scenario numbers in the SAA modelling results are " + negScenarioNum);
			
			// Dispose of model and environment
			model.dispose();
			env.dispose();
		} catch (GRBException e) {
			//System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			result[0] = 0;
			result[1] = 0;
			result[2] = 0;
		}
		return result;
	}
	
	/**
	 * solve the chance SAA by sorting the scenarios,
	 * this is an extended formulation:sort the scenario in each period.
	 * there is a binary variable in each period and each scenario.
	 * bidimap to map the value and key one-one correspondence;
	 * there is delta inequality constraint.
	 * @return
	 */
	public double[] solveSortEach() {
		double[] result = new double[3];
		try {
			int T = distributions.length;
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			int p = negativeScenarioNumRequire; //sampleNumTotal; 
			
			List<List<Integer>> Indexes = new ArrayList<>();
			for (int t = 0; t < T; t++)
				Indexes.add(IntStream.iterate(0, i -> i + 1).limit(sampleNums[t]).boxed().collect(Collectors.toList()));
			CartesianProduct<Integer> CP = new CartesianProduct<>();
			List<List<Integer>> scenarioIndexes = CP.product(Indexes);
		
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);
			//env.set("logFile", "mip-chance.log");
			env.start();

			// Create empty model
			GRBModel model = new GRBModel(env);
			model.set(GRB.IntParam.LogToConsole, 0);

			// Create variables
			GRBVar[][] Q = new GRBVar[T][sampleNumTotal];
			GRBVar[][] delta = new GRBVar[T][sampleNumTotal]; // auxiliary variable
			GRBVar[][] I = new GRBVar[T][sampleNumTotal];
			GRBVar[] beta = new GRBVar[sampleNumTotal];	// whether lost sale happens
			GRBVar[][] alpha = new GRBVar[T][sampleNumTotal];	// auxiliary variable
		    GRBVar[] z = new GRBVar[sampleNumTotal];	// auxiliary variable
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	for (int t = 0; t < T; t++) {
		    		Q[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");		    		
		    		I[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "I");	
		    		alpha[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "alpha");
		    	}
		    	z[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "z");
		    	beta[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "beta");
		    }
		    
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	for (int t = 0; t < T; t++) {
		    		delta[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "delta");
		    	}	
		    }

			// expression variables
			GRBLinExpr[][] cash = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr[][] minCash = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr[][] revenue = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr expectFinalCash = new GRBLinExpr();		
			
			M2 = iniCash + price[0] * M1; // prices are same
			M3 = holdCostUnit * T * iniI + variCostUnit[0] * M1 + Arrays.stream(overheadCost).sum() - iniCash;
			
			// cash flow
		    for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		revenue[t][s] = new GRBLinExpr();
		    		cash[t][s] = new GRBLinExpr();
		    		minCash[t][s] = new GRBLinExpr();
		    		if (t == 0 && T > 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}
		    		else if (t == 0 && T == 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else if (t == T - 1 && T > 1) {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);	
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}	    			
		    	}
		    
		    // objective function, set objective
		    GRBLinExpr expectFinalPositiveCashNum = new GRBLinExpr();
		    expectFinalPositiveCashNum.clear();
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	expectFinalPositiveCashNum.addTerm(1.0, z[s]); // make objective greater 
		    }
			model.setObjective(expectFinalPositiveCashNum, GRB.MAXIMIZE);
			
			// choose maximum sum of demand in T periods as M1
			M1 = 0;
			for (int s = 0; s < sampleNumTotal; s++) {
				double thissSumD = 0;
				for (int t = 0; t < T; t++) {
					int sIndex = scenarioIndexes.get(s).get(t);
					thissSumD += scenarios[t][sIndex];
				}
				if (thissSumD > M1)
					M1 = thissSumD;
			}
						
			// Add constraint
			// inventory flow			
			for (int t = 0; t < T; t++) {	
				
				Integer[] tScenarioIndex = IntStream.iterate(0, i->i+1).limit(sampleNumTotal).boxed().toArray(Integer[]::new);
				Integer[] tScenarioSortIndex = IntStream.iterate(0, i->i+1).limit(sampleNumTotal).boxed().toArray(Integer[]::new);			
				// sort the scenarios in each period, return the index of sorted scenarios
				final int period = t;
				Comparator<Integer> comparator = (o1, o2) -> {
					int thissSumD1 = 0; int thissSumD2 = 0;
					for (int m = 0; m <= period; m++) {
						double demand1 = scenarios[m][scenarioIndexes.get(o1).get(m)];
						thissSumD1 += demand1;
						double demand2 = scenarios[m][scenarioIndexes.get(o2).get(m)];
						thissSumD2 += demand2;
					}
					if (thissSumD1 < thissSumD2) // descending
						return 1;
					else if (thissSumD1 > thissSumD2)// 返回负数表示排在上面，返回正数表示排在下面
						return -1;
					else {
						return 0;
					}
				};				
				Arrays.sort(tScenarioSortIndex, comparator);
				BidiMap<Integer, Integer> bMap = new DualHashBidiMap<Integer, Integer>();
				for (int s = 0; s < sampleNumTotal; s++) {
					int index = tScenarioSortIndex[s];
					bMap.put(index, s);
				}
				
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		int sIndex = scenarioIndexes.get(s).get(t);
		    		double demand = scenarios[t][sIndex];
		    		int thisSSortIndex = bMap.get(s);
		    		if (thisSSortIndex < p) { // delta inequality constraint
		    			if (p < sampleNumTotal) {
		    				int index2 = bMap.getKey(bMap.get(s)+1);		    			
		    				model.addConstr(delta[t][s], GRB.GREATER_EQUAL, delta[t][index2], "deltaConstraint");
		    			}
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
		    		else {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
						if (t==0) {
			    			rightExpr1.addConstant(iniI); rightExpr1.addTerm(1, Q[t][s]); 
			    			rightExpr1.addConstant(-demand); 
						}
						else {
							rightExpr1.addTerm(1, I[t-1][s]); rightExpr1.addTerm(1, Q[t][s]); 
			    			rightExpr1.addConstant(-demand); 	
						}
						model.addConstr(I[t][s], GRB.EQUAL, rightExpr1, "IConstraint1");	
						model.addConstr(delta[t][s], GRB.EQUAL, 0, "deltaEqualZero");
					}
		    	}
			}
			
			// chance constraint
			for (int s = 0; s < sampleNumTotal; s++) {
				for (int t = 0; t < T; t++) {
					model.addConstr(delta[t][s], GRB.LESS_EQUAL, beta[s], "ChanceConstraint1");
				}
			}

			GRBLinExpr sumBeta = new GRBLinExpr();
			for (int s = 0; s < sampleNumTotal; s++) {
				sumBeta.addTerm(1, beta[s]);
			}
			model.addConstr(sumBeta, GRB.LESS_EQUAL, negativeScenarioNumRequire, "ChanceConstraint2");

			// positive cash balance
			M2 = iniCash + price[0] * M1; // prices are same
			M3 = variCostUnit[0] * M1 + Arrays.stream(overheadCost).sum() - iniCash;
						
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++){
//		    		model.addGenConstrIndicator(alpha[t][s], 1, cash[t][s], GRB.GREATER_EQUAL, 0, null);
//		    		model.addGenConstrIndicator(alpha[t][s], 0, cash[t][s], GRB.LESS_EQUAL, 0, null);
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(M2, alpha[t][s]);
		    		model.addConstr(cash[t][s], GRB.LESS_EQUAL, rightExpr, "CashConstraint1");
		    		rightExpr.clear();
		    		rightExpr.addConstant(-M3); rightExpr.addTerm(M3, alpha[t][s]);
		    		model.addConstr(cash[t][s], GRB.GREATER_EQUAL, rightExpr, "CashConstraint4");
		    	}
			
			for (int s = 0; s < sampleNumTotal; s++) {
				GRBLinExpr leftExpr = new GRBLinExpr();
    			leftExpr.addTerm(-1, z[s]); leftExpr.addConstant(1);
    			GRBLinExpr rightExpr = new GRBLinExpr();
    			//rightExpr.clear();
		    	for (int t = 0; t < T; t++){		    
		    		model.addConstr(z[s], GRB.LESS_EQUAL, alpha[t][s], "CashConstraint3");
		    		rightExpr.addTerm(-1, alpha[t][s]);  rightExpr.addConstant(1);
		    	}
		    	model.addConstr(leftExpr, GRB.LESS_EQUAL, rightExpr, "CashConstraint2");
			}
			
			// first-stage decision, here and now decision
			for (int s = 0; s < sampleNumTotal - 1; s++) {
				model.addConstr(Q[0][s], GRB.EQUAL, Q[0][s+1], "fistStageConstraint");
			}
			
			// Optimize model
		    model.optimize();
		    
		    int negScenarioNum = 0;
		    double[] zValue = new double[sampleNumTotal];
		    double[] betaValue = new double[sampleNumTotal];
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
					int thissNegative = 0; 
					for (int t = 0; t < T; t++) {
						out.write(String.format("%.1f  ", cash[t][s].getValue()));
						if (cash[t][s].getValue() < -1)
							thissNegative = 1;
					}
					negScenarioNum += thissNegative;
					zValue[s] = z[s].get(GRB.DoubleAttr.X);
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("alpha in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", alpha[t][s].get(GRB.DoubleAttr.X)));
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("z in each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write(String.format("%.1f  ", z[s].get(GRB.DoubleAttr.X)));	
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
			    	betaValue[s] = beta[s].get(GRB.DoubleAttr.X);
			    }
				out.close();
			} catch (IOException e) {				
			}	    
			
			result[0] = Q[0][0].get(GRB.DoubleAttr.X);
			result[1] = model.get(GRB.DoubleAttr.ObjVal);
			result[2] = Arrays.stream(betaValue).sum();
			//System.out.println("negative cash scenario numbers in the SAA modelling results are " + negScenarioNum);
			
			// Dispose of model and environment
			model.dispose();
			env.dispose();
		} catch (GRBException e) {
			//System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			result[0] = 0;
			result[1] = 0;
			result[2] = 0;
		}
		return result;
	}
	
	
	
	/**
	 * sort the scenarios in the whole planning horizon, not in the single period.
	 * @return
	 */
	public double[] solveSortWhole() {
		double[] result = new double[3];
		try {
			int T = distributions.length;
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			int p = negativeScenarioNumRequire; //sampleNumTotal; 
			
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
				if (sum1 < sum2) // 返回负数表示排在上面，返回正数表示排在下面
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
			//env.set("logFile", "mip-chance.log");
			env.start();

			// Create empty model
			GRBModel model = new GRBModel(env);
			model.set(GRB.IntParam.LogToConsole, 0);

			// Create variables
			GRBVar[][] Q = new GRBVar[T][sampleNumTotal];
			GRBVar[][] delta = new GRBVar[T][p]; // auxiliary variable, only p
			GRBVar[][] I = new GRBVar[T][sampleNumTotal];
			GRBVar[] beta = new GRBVar[sampleNumTotal];	// whether lost sale happens
			GRBVar[][] alpha = new GRBVar[T][sampleNumTotal];	// auxiliary variable
		    GRBVar[] z = new GRBVar[sampleNumTotal];	// auxiliary variable
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	for (int t = 0; t < T; t++) {
		    		Q[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");		    		
		    		I[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "I");	
		    		alpha[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "alpha");
		    	}
		    	z[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "z");
		    	beta[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "beta");
		    }
		    
		    for (int t = 0; t < T; t++) {
		    	for (int s = 0; s < p; s++) {
		    		delta[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "delta");	
		    	}
		    }

			// expression variables
			GRBLinExpr[][] cash = new GRBLinExpr[T][sampleNumTotal];
//			GRBLinExpr[][] minCash = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr[][] revenue = new GRBLinExpr[T][sampleNumTotal];
//			GRBLinExpr expectFinalCash = new GRBLinExpr();		
			
			// choose maximum sum of demand in T periods as M1
			M1 = 0;
			try {
			for (int s = 0; s < sampleNumTotal; s++) {
				double thissSumD = 0;
				for (int t = 0; t < T; t++) {
					thissSumD += listScenariosAll.get(s).get(t);
				}
				if (thissSumD > M1)
					M1 = thissSumD;
			
			}
				
			} catch (Exception e) {
				System.out.println();
			}
			
			
			M2 = iniCash + price[0] * M1; // prices are same
			M3 = holdCostUnit * T * iniI + variCostUnit[0] * M1 + Arrays.stream(overheadCost).sum() - iniCash;
			
			// cash flow
		    for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		revenue[t][s] = new GRBLinExpr();
		    		cash[t][s] = new GRBLinExpr();
//		    		minCash[t][s] = new GRBLinExpr();
//		    		if (s < p) {
			    		if (t == 0 && T > 1) {
			    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
			    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
			    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
			    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
			    		}
			    		else if (t == 0 && T == 1) {
			    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
			    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
			    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
			    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
			    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
						}
			    		else if (t == T - 1 && T > 1) {
			    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
			    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
			    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
			    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
			    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
						}
			    		else {
			    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
			    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
			    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);	
			    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
			    		}	    	
		    		}
//		    		else { // 库存持有成本计算时还需要修改，要累加
//		    			GRBLinExpr HIts = new GRBLinExpr();
//		    			GRBLinExpr salValue = new GRBLinExpr();
//		    			cash[t][s].addConstant(iniCash); 
////		    			HIts.addConstant(-holdCostUnit*iniI); 
//		    			salValue.addConstant(salvageValueUnit*iniI); 
//		    			double pD = 0;
//		    			for (int k = 0; k <= t; k++) {
//		    				double demand = listScenariosAll.get(s).get(k);
//		    				pD = price[k]*demand;
//		    				cash[t][s].addConstant(pD);
//		    				cash[t][s].addTerm(-variCostUnit[k], Q[k][s]);
//		    				cash[t][s].addConstant(-overheadCost[k]);
////		    				HIts.addTerm(-holdCostUnit, Q[k][s]);
////		    				HIts.addConstant(demand*holdCostUnit);
//		    				salValue.addTerm(salvageValueUnit, Q[k][s]);
//		    				salValue.addConstant(-demand*salvageValueUnit);
//		    			}
////		    			cash[t][s].add(HIts);
//		    			if (t == T-1)
//		    				cash[t][s].add(salValue);
//					}
//		    	}
		    
		    // objective function, set objective
		    GRBLinExpr expectFinalPositiveCashNum = new GRBLinExpr();
		    expectFinalPositiveCashNum.clear();
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	expectFinalPositiveCashNum.addTerm(1.0, z[s]); // make objective greater 
		    }
			model.setObjective(expectFinalPositiveCashNum, GRB.MAXIMIZE);
			

						
			// Add constraint
			// inventory flow			
			for (int t = 0; t < T; t++) {					
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		double demand = listScenariosAll.get(s).get(t);
		    		if (s < p) {
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
		    		else {// 服务率高，下面这些等式约束条件就多，反而总体计算慢了
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
						if (t==0) {
			    			rightExpr1.addConstant(iniI); rightExpr1.addTerm(1, Q[t][s]); 
			    			rightExpr1.addConstant(-demand); 
						}
						else {
							rightExpr1.addTerm(1, I[t-1][s]); rightExpr1.addTerm(1, Q[t][s]); 
			    			rightExpr1.addConstant(-demand); 	
						}
						model.addConstr(I[t][s], GRB.EQUAL, rightExpr1, "IConstraint1");	
					}
		    	}
			}
			
			// chance constraint
			for (int s = 0; s < p; s++) {
				for (int t = 0; t < T; t++) {
					model.addConstr(delta[t][s], GRB.LESS_EQUAL, beta[s], "ChanceConstraint1");
				}
			}

			GRBLinExpr sumBeta = new GRBLinExpr();
			for (int s = 0; s < sampleNumTotal; s++) {
				sumBeta.addTerm(1, beta[s]);
			}
			model.addConstr(sumBeta, GRB.LESS_EQUAL, negativeScenarioNumRequire, "ChanceConstraint2");

			// positive cash balance
			M2 = iniCash + price[0] * M1; // prices are same
			M3 = variCostUnit[0] * M1 + Arrays.stream(overheadCost).sum() - iniCash;
						
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++){
//		    		model.addGenConstrIndicator(alpha[t][s], 1, cash[t][s], GRB.GREATER_EQUAL, 0, null);
//		    		model.addGenConstrIndicator(alpha[t][s], 0, cash[t][s], GRB.LESS_EQUAL, 0, null);
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(M2, alpha[t][s]);
		    		model.addConstr(cash[t][s], GRB.LESS_EQUAL, rightExpr, "CashConstraint1");
		    		rightExpr.clear();
		    		rightExpr.addConstant(-M3); rightExpr.addTerm(M3, alpha[t][s]);
		    		model.addConstr(cash[t][s], GRB.GREATER_EQUAL, rightExpr, "CashConstraint4");
		    	}
			
			for (int s = 0; s < sampleNumTotal; s++) {
				GRBLinExpr leftExpr = new GRBLinExpr();
    			leftExpr.addTerm(-1, z[s]); leftExpr.addConstant(1);
    			GRBLinExpr rightExpr = new GRBLinExpr();
    			//rightExpr.clear();
		    	for (int t = 0; t < T; t++){		    
		    		model.addConstr(z[s], GRB.LESS_EQUAL, alpha[t][s], "CashConstraint3");
		    		rightExpr.addTerm(-1, alpha[t][s]);  rightExpr.addConstant(1);
		    	}
		    	model.addConstr(leftExpr, GRB.LESS_EQUAL, rightExpr, "CashConstraint2");
			}
			
			// first-stage decision, here and now decision
			for (int s = 0; s < sampleNumTotal - 1; s++) {
				model.addConstr(Q[0][s], GRB.EQUAL, Q[0][s+1], "fistStageConstraint");
			}
			
			// Optimize model
		    model.optimize();
		    
		    int negScenarioNum = 0;
		    double[] zValue = new double[sampleNumTotal];
		    double[] betaValue = new double[sampleNumTotal];
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
					int thissNegative = 0; 
					for (int t = 0; t < T; t++) {
						out.write(String.format("%.1f  ", cash[t][s].getValue()));
						if (cash[t][s].getValue() < -1)
							thissNegative = 1;
					}
					negScenarioNum += thissNegative;
					zValue[s] = z[s].get(GRB.DoubleAttr.X);
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("alpha in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", alpha[t][s].get(GRB.DoubleAttr.X)));
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("z in each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write(String.format("%.1f  ", z[s].get(GRB.DoubleAttr.X)));	
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
			    	betaValue[s] = beta[s].get(GRB.DoubleAttr.X);
			    }
				out.close();
			} catch (IOException e) {				
			}	    
			
			result[0] = Q[0][0].get(GRB.DoubleAttr.X);
			result[1] = model.get(GRB.DoubleAttr.ObjVal);
			result[2] = Arrays.stream(betaValue).sum();
			//System.out.println("negative cash scenario numbers in the SAA modelling results are " + negScenarioNum);
			
			// Dispose of model and environment
			model.dispose();
			env.dispose();
		} catch (GRBException e) {
			//System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			result[0] = 0;
			result[1] = 0;
			result[2] = 0;
		}
		return result;
	}
	
	
	/**
	 * solve the chance SAA by sorting the scenarios in each period,
	 * this is an extended formulation.
	 * sort the scenario in each period but without the delta inequality constraint.
	 * @return
	 */
	public double[] solveSort2() {
		double[] result = new double[3];
		try {
			int T = distributions.length;
			int sampleNumTotal = IntStream.of(sampleNums).reduce(1, (a, b) -> a * b);
			int negativeScenarioNumRequire = (int) (sampleNumTotal * (1 - serviceRate));
			int p = negativeScenarioNumRequire; //sampleNumTotal; 
			
	
			// use Gurobi to solve the mip model
			// Create empty environment, set options, and start
			GRBEnv env = new GRBEnv(true);
			//env.set("logFile", "mip-chance.log");
			env.start();

			// Create empty model
			GRBModel model = new GRBModel(env);
			model.set(GRB.IntParam.LogToConsole, 0);

			// Create variables
			GRBVar[][] Q = new GRBVar[T][sampleNumTotal];
			GRBVar[][] delta = new GRBVar[T][p]; // auxiliary variable, only p
			GRBVar[][] I = new GRBVar[T][sampleNumTotal];
			GRBVar[] beta = new GRBVar[sampleNumTotal];	// whether lost sale happens
			GRBVar[][] alpha = new GRBVar[T][sampleNumTotal];	// auxiliary variable
		    GRBVar[] z = new GRBVar[sampleNumTotal];	// auxiliary variable
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	for (int t = 0; t < T; t++) {
		    		Q[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "Q");		    		
		    		I[t][s] = model.addVar(0.0, Double.MAX_VALUE, 0.0, GRB.CONTINUOUS, "I");	
		    		alpha[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "alpha");
		    	}
		    	z[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "z");
		    	beta[s] = model.addVar(0.0, 1, 0.0, GRB.BINARY, "beta");
		    }
		    
		    for (int t = 0; t < T; t++) {
		    	for (int s = 0; s < p; s++) {
		    		delta[t][s] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "delta");	
		    	}
		    }

			// expression variables
			GRBLinExpr[][] cash = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr[][] minCash = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr[][] revenue = new GRBLinExpr[T][sampleNumTotal];
			GRBLinExpr expectFinalCash = new GRBLinExpr();		
			
			M2 = iniCash + price[0] * M1; // prices are same
			M3 = holdCostUnit * T * iniI + variCostUnit[0] * M1 + Arrays.stream(overheadCost).sum() - iniCash;
			
			// cash flow
		    for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		revenue[t][s] = new GRBLinExpr();
		    		cash[t][s] = new GRBLinExpr();
		    		minCash[t][s] = new GRBLinExpr();
		    		if (t == 0 && T > 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}
		    		else if (t == 0 && T == 1) {
		    			revenue[t][s].addConstant(price[t] * iniI); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].addConstant(iniCash); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else if (t == T - 1 && T > 1) {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);
		    			cash[t][s].addTerm(salvageValueUnit, I[t][s]);
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
					}
		    		else {
		    			revenue[t][s].addTerm(price[t], I[t-1][s]); revenue[t][s].addTerm(price[t], Q[t][s]); revenue[t][s].addTerm(-price[t], I[t][s]);
		    			cash[t][s].multAdd(1, cash[t-1][s]); cash[t][s].add(revenue[t][s]); 
		    			cash[t][s].addTerm(-variCostUnit[t], Q[t][s]); cash[t][s].addConstant(-overheadCost[t]);	
		    			cash[t][s].addTerm(-holdCostUnit, I[t][s]);
		    		}	    			
		    	}
		    
		    // objective function, set objective
		    GRBLinExpr expectFinalPositiveCashNum = new GRBLinExpr();
		    expectFinalPositiveCashNum.clear();
		    for (int s = 0; s < sampleNumTotal; s++) {
		    	expectFinalPositiveCashNum.addTerm(1.0, z[s]); // make objective greater 
		    }
			model.setObjective(expectFinalPositiveCashNum, GRB.MAXIMIZE);
			
			// choose maximum sum of demand in T periods as M1
			M1 = 10000;
	
						
			// Add constraint
			// inventory flow		
			// sort scenario directly in each period, does not use bmap to store index of the sorted scenario
			for (int t = 0; t < T; t++) {		
				
				List<List<Double>> listScenarios = new ArrayList<>();
				for (int t1 = 0; t1 < t + 1; t1++)
					listScenarios.add(Arrays.stream(scenarios[t1]).boxed().collect(Collectors.toList()));
				CartesianProduct<Double> CP = new CartesianProduct<>();
				List<List<Double>> listScenariosAll = CP.product(listScenarios);
				Comparator<List<Double>> comparator = (o1, o2) -> {
					double sum1 = 0;
					double sum2 = 0;
					for (int t1 = 0; t1 < o1.size(); t1++) {
						sum1 += o1.get(t1)*price[t1];
						sum2 += o2.get(t1)*price[t1];
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
				
				int divideNum = IntStream.of(sampleNums).skip(t + 1).reduce(1, (x, y) -> x*y);
		    	for (int s = 0; s < sampleNumTotal; s++) {
		    		// add the s index
		    		int sIndex = s / divideNum;
		    		double demand = listScenariosAll.get(sIndex).get(t);
		    		if (s < p) {
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
		    		else {
		    			GRBLinExpr rightExpr1 = new GRBLinExpr();
						if (t==0) {
			    			rightExpr1.addConstant(iniI); rightExpr1.addTerm(1, Q[t][s]); 
			    			rightExpr1.addConstant(-demand); 
						}
						else {
							rightExpr1.addTerm(1, I[t-1][s]); rightExpr1.addTerm(1, Q[t][s]); 
			    			rightExpr1.addConstant(-demand); 	
						}
						model.addConstr(I[t][s], GRB.EQUAL, rightExpr1, "IConstraint1");	
//						model.addConstr(delta[t][s], GRB.EQUAL, 0, "deltaEqualZero");
					}
		    	}
			}
			
			// chance constraint
			for (int s = 0; s < p; s++) {
				for (int t = 0; t < T; t++) {
					model.addConstr(delta[t][s], GRB.LESS_EQUAL, beta[s], "ChanceConstraint1");
				}
			}

			GRBLinExpr sumBeta = new GRBLinExpr();
			for (int s = 0; s < sampleNumTotal; s++) {
				sumBeta.addTerm(1, beta[s]);
			}
			model.addConstr(sumBeta, GRB.LESS_EQUAL, negativeScenarioNumRequire, "ChanceConstraint2");

			// positive cash balance
			M2 = iniCash + price[0] * M1; // prices are same
			M3 = variCostUnit[0] * M1 + Arrays.stream(overheadCost).sum() - iniCash;
						
			for (int t = 0; t < T; t++) 
		    	for (int s = 0; s < sampleNumTotal; s++){
//		    		model.addGenConstrIndicator(alpha[t][s], 1, cash[t][s], GRB.GREATER_EQUAL, 0, null);
//		    		model.addGenConstrIndicator(alpha[t][s], 0, cash[t][s], GRB.LESS_EQUAL, 0, null);
		    		GRBLinExpr rightExpr = new GRBLinExpr();
		    		rightExpr.addTerm(M2, alpha[t][s]);
		    		model.addConstr(cash[t][s], GRB.LESS_EQUAL, rightExpr, "CashConstraint1");
		    		rightExpr.clear();
		    		rightExpr.addConstant(-M3); rightExpr.addTerm(M3, alpha[t][s]);
		    		model.addConstr(cash[t][s], GRB.GREATER_EQUAL, rightExpr, "CashConstraint4");
		    	}
			
			for (int s = 0; s < sampleNumTotal; s++) {
				GRBLinExpr leftExpr = new GRBLinExpr();
    			leftExpr.addTerm(-1, z[s]); leftExpr.addConstant(1);
    			GRBLinExpr rightExpr = new GRBLinExpr();
    			//rightExpr.clear();
		    	for (int t = 0; t < T; t++){		    
		    		model.addConstr(z[s], GRB.LESS_EQUAL, alpha[t][s], "CashConstraint3");
		    		rightExpr.addTerm(-1, alpha[t][s]);  rightExpr.addConstant(1);
		    	}
		    	model.addConstr(leftExpr, GRB.LESS_EQUAL, rightExpr, "CashConstraint2");
			}
			
			// first-stage decision, here and now decision
			for (int s = 0; s < sampleNumTotal - 1; s++) {
				model.addConstr(Q[0][s], GRB.EQUAL, Q[0][s+1], "fistStageConstraint");
			}
			
			// Optimize model
		    model.optimize();
		    
		    int negScenarioNum = 0;
		    double[] zValue = new double[sampleNumTotal];
		    double[] betaValue = new double[sampleNumTotal];
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
					int thissNegative = 0; 
					for (int t = 0; t < T; t++) {
						out.write(String.format("%.1f  ", cash[t][s].getValue()));
						if (cash[t][s].getValue() < -1)
							thissNegative = 1;
					}
					negScenarioNum += thissNegative;
					zValue[s] = z[s].get(GRB.DoubleAttr.X);
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("alpha in each period and each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write("scenario " + s + ": \n");
					for (int t = 0; t < T; t++) 
						out.write(String.format("%.1f  ", alpha[t][s].get(GRB.DoubleAttr.X)));
					out.newLine();
				}
				out.newLine();
				out.newLine();
				out.newLine();
				out.write("***************************************************\n");
				out.write("z in each scenario: \n");
				for (int s = 0; s < sampleNumTotal; s++) {
					out.write(String.format("%.1f  ", z[s].get(GRB.DoubleAttr.X)));	
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
			    	betaValue[s] = beta[s].get(GRB.DoubleAttr.X);
			    }
				out.close();
			} catch (IOException e) {				
			}	    
			
			result[0] = Q[0][0].get(GRB.DoubleAttr.X);
			result[1] = model.get(GRB.DoubleAttr.ObjVal);
			result[2] = Arrays.stream(betaValue).sum();
			//System.out.println("negative cash scenario numbers in the SAA modelling results are " + negScenarioNum);
			
			// Dispose of model and environment
			model.dispose();
			env.dispose();
		} catch (GRBException e) {
			//System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
			result[0] = 0;
			result[1] = 0;
			result[2] = 0;
		}
		return result;
	}
	
}
	
	