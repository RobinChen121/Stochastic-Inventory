package milp;

import java.util.Arrays;

import milp.BinaryMILP.BoundCriteria;
import sdp.write.WriteToCsv;

/**
* @author Zhen Chen
* @date: 2018年11月14日 下午7:55:11  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  testing the performance of binary MILP
* 
* my testing shows its average gap with sdp is 5%
* after adding extra piecewise constraints, the average gap is 
*/

public class BinaryMILPTesting {

	public static void main(String[] args) {
		String headString = "K, v, h, I0, pai, coeVar, DemandPatt, simMILPValue, Time(sec)";
		WriteToCsv.writeToFile("./" + "testMILP_results.csv", headString);
		
		double[][] demands =
			 {{10,	10,	10,	10,	10,	10,	10,	10},
			 {15,	16,	15,	14,	11,	7,	6,	3},
			 {3,	6,	7,	11,	14,	15,	16,	15},
			 {15,   4,	4,	10,	18,	 4,	 4,	10},
			 {12,	7,	7,	10,	13,	7,	7,	12},
			 {2,	4,	7,	3,	10,	10,	3,	3},
			 {5,	15,	26,	44,	24,	15,	22,	10},
			 {4,	23,	28,	50,	39,	26,	19,	32},
			 {11,	14,	7,	11,	16,	31,	11,	48},
			 {18,	6,	22,	22,	51,	54,	22,	21},
		};
		
		double[] K = { 200, 300, 400 }; 
		double[] v = { 0, 1 };
		double[] pai = { 5, 10, 20};		
		double[] coeVar = {0.1, 0.2, 0.3};
		double holdingCost = 1;
		double iniInventory = 0;
		
		for (int idemand = 3; idemand < 10; idemand++)
			for (int iv = 0; iv < v.length; iv++) 
				for (int ipai = 0; ipai < pai.length; ipai++) 
					for (int iK = 0; iK < K.length; iK++) 
						for (int iCoe = 0; iCoe < coeVar.length; iCoe++) {
							double[] meanDemand = demands[idemand];
							double fixOrderCost = K[iK];
							double variCost = v[iv];
							double penaltyCost = pai[ipai];
							double coeVarValue = coeVar[iCoe];
							double[] sigma = Arrays.stream(meanDemand).map(i -> coeVarValue*i).toArray();
							int partionNum = 10;
							BoundCriteria boundCriteria = BoundCriteria.LOWBOUND;
							
							long currTime = System.currentTimeMillis();
							BinaryMILP binary = new BinaryMILP(meanDemand, sigma, iniInventory, fixOrderCost, variCost, holdingCost, penaltyCost, partionNum, boundCriteria);
							double[][] sS = binary.binaryFind();
							System.out.println("Approximate sS are: ");
							System.out.println(Arrays.deepToString(sS));
							double time = (System.currentTimeMillis() - currTime) / 1000;
							System.out.println("running time is " + time + "s");
							
							/*******************************************************************
							 * Simulating BinaryApprox method
							 */
							int sampleNum = 10000;
							Simulation simulate = new Simulation(meanDemand, sigma, iniInventory, fixOrderCost, variCost, holdingCost, penaltyCost);
							double simMILPlValue = simulate.simulatesS(sS, sampleNum);
							System.out.println("***************************************************");
							
							String out = fixOrderCost + ",\t" 
									+ variCost + ",\t" 
									+ holdingCost + ",\t" 
									+ iniInventory + ",\t" 
									+ penaltyCost + ",\t" 
									+ coeVar[iCoe] + ",\t"
									+ (idemand + 1) + ",\t" 
									+ simMILPlValue + ",\t" 
									+ time;
							WriteToCsv.writeToFile("./" + "testMILP_results.csv", out);						
						}
	}
}


