package capacitated;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import sdp.inventory.GetPmf;
import sdp.inventory.State;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 9, 2018---12:48:56 PM
 * @description: solving capacitated lot sizing problem, the code imitate
 *               Roberto Rossi's; use Interface to define several functions and
 *               apply lambda expressions
 * 
 */

public class CLSP {
	double[][][] pmf;

	public CLSP(double[][][] pmf) {
		this.pmf = pmf;
	}

	Function<State, double[]> actionGenerator;

	@FunctionalInterface
	interface StateTransitionFunction<S, A, R, S2> {
		public S2 apply(S s, A a, R r);
	}

	StateTransitionFunction<State, Double, Double, State> stateTransition;

	@FunctionalInterface
	interface ImmediateValueFunction<S, A, R, V> {
		public V apply(S s, A a, R r);
	}

	ImmediateValueFunction<State, Double, Double, Double> immediateValue;

	Map<State, Double> cacheValues = new HashMap<>();
	SortedMap<State, Double> cacheActions = new TreeMap<>();

	double f(State state) {
		return cacheValues.computeIfAbsent(state, s -> {
			double val = Arrays.stream(actionGenerator.apply(s)).map(orderQty -> Arrays.stream(pmf[s.getPeriod() - 1])
					.parallel()
					.mapToDouble(p -> p[1] * immediateValue.apply(s, orderQty, p[0])
							+ (s.getPeriod() < pmf.length ? p[1] * f(stateTransition.apply(s, orderQty, p[0])) : 0))
					.sum()).min().getAsDouble();
			double bestOrderQty = Arrays.stream(actionGenerator.apply(s)).filter(orderQty -> Arrays
					.stream(pmf[s.getPeriod() - 1])
					.mapToDouble(p -> p[1] * immediateValue.apply(s, orderQty, p[0])
							+ (s.getPeriod() < pmf.length ? p[1] * f(stateTransition.apply(s, orderQty, p[0])) : 0))
					.sum() == val).findAny().getAsDouble();
			cacheActions.putIfAbsent(s, bestOrderQty);
			return val;
		});
	}

	int[] levelNum(double[][] optTable, int maxOrderQuantity) {
		ArrayList<Integer> indexArr = new ArrayList<>();
		boolean mark = false;
		for (int j = 0; j < optTable.length; j++) {
			if (optTable[j][2] < maxOrderQuantity && !mark) {
				mark = true;
			} else if (optTable[j][2] == maxOrderQuantity && mark) {
				mark = false;
				indexArr.add(j);
			}
			if (optTable[j][2] == 0) {
				indexArr.add(j);
				break;
			}
		}
		return indexArr.stream().mapToInt(p -> p.intValue()).toArray();
	}

	public static void main(String[] args) {
		double initialInventory = 0;
		double[] meanDemand = { 9, 23, 53, 29 };

		double truncationQuantile = 0.9999;
		double stepSize = 1;
		double minState = -200;
		double maxState = 300;
		int maxOrderQuantity = 100;

		double fixedOrderingCost = 500;
		double proportionalOrderingCost = 0;
		double penaltyCost = 10;
		double holdingCost = 2;

		int T = meanDemand.length;
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(T)
				.mapToObj(i -> new PoissonDist(meanDemand[i])).toArray(PoissonDist[]::new);
		double[][][] pmf = new GetPmf(distributions, truncationQuantile, stepSize).getpmf();

		CLSP inventory = new CLSP(pmf);

		inventory.actionGenerator = s -> {
			return DoubleStream.iterate(0, i -> i + stepSize).limit(maxOrderQuantity + 1).toArray();
		};

		inventory.stateTransition = (state, action, randomDemand) -> {
			double nextInventory = state.getIniInventory() + action - randomDemand;
			nextInventory = nextInventory > maxState ? maxState : nextInventory;
			nextInventory = nextInventory < minState ? minState : nextInventory;
			return new State(state.getPeriod() + 1, nextInventory);
		};

		inventory.immediateValue = (state, action, randomDemand) -> {
			double fixedCost = action > 0 ? fixedOrderingCost : 0;
			double variableCost = proportionalOrderingCost * action;
			double inventoryLevel = state.getIniInventory() + action - randomDemand;
			double holdingCosts = holdingCost * Math.max(inventoryLevel, 0);
			double penaltyCosts = penaltyCost * Math.max(-inventoryLevel, 0);
			double totalCosts = fixedCost + variableCost + holdingCosts + penaltyCosts;
			return totalCosts;
		};

		int period = 1;
		State initialState = new State(period, initialInventory);
		long currTime2 = System.currentTimeMillis();

		double finalValue = inventory.f(initialState);
		System.out.println("final optimal expected value is: " + finalValue);
		double optQ = inventory.cacheActions.get(new State(period, initialInventory));
		System.out.println("optimal order quantity in the first priod is : " + optQ);
		double time = (System.currentTimeMillis() - currTime2) / 1000;
		System.out.println("running time is " + time + "s");
	}
}
