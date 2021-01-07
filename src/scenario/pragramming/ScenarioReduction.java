/**
 * @date: May 17, 2020
 */
package scenario.pragramming;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;


/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: May 17, 2020
 * @Desc: implement the scenario reduction method in Hu and Hu (2016)
 * 
 * distances matrix should not be too large, out of memory, and this is same for hashmap
 * 
 * demand follow non-stationary gamma distribution.
               3-item, 6 periods:
                Distribution   & Gamma & Gamma & Gamma  \\
                Scale   & 62.99 & 199.34&124.05\\
                Shape &1.30 & 1.99 &1.46 \\
                Mean &82.14 & 397.33 &181.54 \\
                Variance &5173.98 &79206.22 &22520.71 \\
                Skewness &2.06 &0.41 &1.67 \\
                Kurtosis &7.39 & 1.78 &5.37\\
               
        scenario tree 1: possibility and demand realizations
                    0.105  & 25  & 290  & 109  \\ 
                    0.341  & 58  & 365  & 90  \\ 
                    0.330  & 62  & 134  & 132  \\ 
                    0.106  & 289  & 789  & 273  \\ 
                    0.119  & 74  & 965  & 564  \\
 *
 */
public class ScenarioReduction {

	
	static double euclDistance(double[][] arr, double[] s1, double[] s2) {
		double squreSum = 0;
		for (int k = 0; k < s1.length - 1; k++)
			for (int i = 0; i < 3; i++) {
				int index1 = (int) s1[k];
				int index2 = (int) s2[k];
				squreSum += Math.pow(arr[index1][i] - arr[index2][i], 2);
			}
		return Math.sqrt(squreSum);
	}
	
	
	public static void main(String[] args) {
		//int itemNum = 3;
		
		int K = 15; // K is the final scenario number needed
		int T = 6;
		int realizationNum = 5;
		double[] possiblities = {0.105, 0.341, 0.330, 0.106, 0.119};
		double[][] demandRealizations = {{25, 290, 109}, {58, 365, 90}, {62, 134, 132}, {289, 789, 273}, {74, 965, 564}};
		int scenarioNum = (int) Math.pow(realizationNum, T);
		double[][] scenarioIndexP = new double[scenarioNum][T + 1];
		
		// all the possible scenarios
		int index = 0;
		for (int i1 = 0; i1 < realizationNum; i1++)
			for (int i2 = 0; i2 < realizationNum; i2++)
				for (int i3 = 0; i3 < realizationNum; i3++)
					for (int i4 = 0; i4 < realizationNum; i4++)
						for (int i5 = 0; i5 < realizationNum; i5++)
							for (int i6 = 0; i6 < realizationNum; i6++) {
								double ratio = possiblities[i1]*possiblities[i2]*possiblities[i3]*possiblities[i4]*possiblities[i5]*possiblities[i6];
								scenarioIndexP[index] = new double[]{i1, i2, i3, i4, i5, i6, ratio};
								index++;
							}

		// select K scenarios
		int selected = 0;
		ArrayList<Integer> indexRecord = new ArrayList<>();
		ArrayList<Double> possibRecord = new ArrayList<>();
		while (selected < K) {
			double[] weightedDistance = new double[scenarioNum];
			double minDistance = 100000;
			int minDIndex = 0;	
			index = 0;
			for (int i = 0; i < scenarioNum; i++) {
				if (!indexRecord.contains(i)) {
					double upperWeightD = 0;
					double lowerWeightD = 0;
					double[] scenarioIndexP1 = new double[T + 1];
					double[] scenarioIndexP2 = new double[T + 1];
					for (int j1 = 0; j1 < i; j1 ++) {
						if (!indexRecord.contains(i)) {
							scenarioIndexP1 = scenarioIndexP[j1]; 
							scenarioIndexP2 = scenarioIndexP[i]; 
							double thisDistance = euclDistance(demandRealizations, scenarioIndexP1, scenarioIndexP2);
							for (int k = 0; k < indexRecord.size(); k++) {
								double[] scenarioIndexPL = scenarioIndexP[indexRecord.get(k)];
								double tempDistance = euclDistance(demandRealizations, scenarioIndexP1, scenarioIndexPL);
								if (tempDistance < thisDistance) 
									thisDistance = tempDistance;
							}
							upperWeightD += scenarioIndexP[j1][6] * thisDistance;
						}
					}
					for (int j2 = i + 1; j2 < scenarioNum; j2++) {
						if (!indexRecord.contains(j2)) {
							scenarioIndexP1 = scenarioIndexP[i]; 
							scenarioIndexP2 = scenarioIndexP[j2]; 
							double thisDistance = euclDistance(demandRealizations, scenarioIndexP1, scenarioIndexP2);
							for (int k = 0; k < indexRecord.size(); k++) {
								double[] scenarioIndexPL = scenarioIndexP[indexRecord.get(k)];
								double tempDistance = euclDistance(demandRealizations, scenarioIndexP2, scenarioIndexPL);
								if (tempDistance < thisDistance) 
									thisDistance = tempDistance;
							}
							lowerWeightD += scenarioIndexP[j2][6] * thisDistance;	
						}
					}
					weightedDistance[i] = upperWeightD + lowerWeightD;
					if (weightedDistance[i] < minDistance) {
						minDistance = weightedDistance[i];
						minDIndex = i;			
					}
					index++;
				}
				else {
					continue;
				}
			}
			indexRecord.add(minDIndex);
			possibRecord.add(scenarioIndexP[minDIndex][T]);
			System.out.println(minDistance);
			selected++;		
		}		
		
		// add possibilities
		for (int i = 0; i < scenarioNum; i++) {
			if (!indexRecord.contains(i)) {
				double[] scenarioIndexP1 = new double[T + 1];
				double[] scenarioIndexP2 = new double[T + 1];
				double minDistance = 100000;
				int minDIndex = 0;
				for (int j = 0; j < indexRecord.size(); j++) {
					scenarioIndexP1 = scenarioIndexP[i]; 
					scenarioIndexP2 = scenarioIndexP[indexRecord.get(j)];
					double thisDistance = euclDistance(demandRealizations, scenarioIndexP1, scenarioIndexP2);
					if (thisDistance < minDistance) {
						minDistance = thisDistance;
						minDIndex = j;
					}
				}
				double newP = possibRecord.get(minDIndex) + scenarioIndexP[i][T];
				possibRecord.set(minDIndex, newP);		

			}
		}
		
		System.out.println(indexRecord.toString());
		System.out.println();
		System.out.println(possibRecord.toString());
		
		// output scenario index
		System.out.println("the scenarios are: ");
		for (int j = 0; j < indexRecord.size(); j++) {
			double[] scenario = scenarioIndexP[indexRecord.get(j)]; 
			System.out.println(Arrays.toString(scenario));
		}
		
	}
	
	
	

}
