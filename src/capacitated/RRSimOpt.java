package capacitated;

import java.util.Arrays;
import java.util.stream.IntStream;

import sdp.sampling.Sampling;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.stat.TallyStore;


/**
 * @author: Roberto Rossi
 * @email: 15011074486@163.com
 * @date: Jul 9, 2018---12:48:56 PM
 * @description: Roberto Rossi's code for Simulation heuristic
 * 
 */


public class RRSimOpt {

	int minMonteCarloSimulationRunsQ = 100;
	int minMonteCarloSimulationRunsP = 100;

	double minQ = 0;
	double maxQ;

	double confidenceLevel = 0.90;
	double percentageError = 1;
	double step = 0.1;

	Sampling sfQ = new Sampling();
	Sampling sfP = new Sampling();

	double K;
	double v;
	double s;
	double h;
	double pai;

	double inventory;

	Distribution[] distributions;

	public RRSimOpt(double K, double v, double pai, double h, Distribution[] distributions, double inventory,
			double maxOrderQuantity) {
		this.K = K;
		this.v = v;
		this.h = h;
		this.pai = pai;
		this.maxQ = maxOrderQuantity;

		this.distributions = distributions;

		this.inventory = inventory;

		Sampling.resetStartStream();
	}

	public double orderingCost(double Q) {
		return Q > 0 ? this.K + Q * this.v : 0;
	}

	public double inventoryCost(double inventory) {
		return this.h * Math.max(inventory, 0) + this.pai * Math.max(-inventory, 0);
	}

	public double simulateSingleRunCycle(double Q, double inventory, Distribution[] distributions) {
		double i = inventory;
		double c = 0;
		double[][] d = Sampling.generateLHSamples(distributions, 1);
		for (int t = 0; t < distributions.length; t++) {
			if (t == 0) {
				c += orderingCost(Q);
				i += Q;
			}

			i -= d[t][0];
			c += inventoryCost(i);
		}
		return c;
	}

	public double[] simulateCycle(double Q, double inventory, Distribution[] distributions) {
		Sampling.resetStartStream();
		TallyStore observationsTally = new TallyStore();
		for (int runs = 0; runs < this.minMonteCarloSimulationRunsQ; runs++) {
			observationsTally.add(simulateSingleRunCycle(Q, inventory, distributions));
		}
		double[] centerAndRadius = new double[2];
		observationsTally.confidenceIntervalStudent(this.confidenceLevel, centerAndRadius);
		/*
		 * while(Math.abs(centerAndRadius[1]/centerAndRadius[0]) >=
		 * percentageError/100){ observationsTally.add(simulateSingleRunCycle(Q,
		 * inventory, capital, distributions));
		 * observationsTally.confidenceIntervalStudent(confidenceLevel,
		 * centerAndRadius); }
		 */
		return centerAndRadius;
	}

	public double getQ(int period, double inventory) {
		double noOrder = simulateCycle(0, inventory, new Distribution[] { this.distributions[period] })[0];

		double bestCapital = Double.MAX_VALUE;
		double bestQ = 0;
		for (int t = 0; t < this.distributions.length - period; t++) {
			Distribution[] reducedHorizon = Arrays.stream(distributions, period, period + t + 1)
					.toArray(Distribution[]::new);

			double Qlb = minQ;
			double Qub = maxQ;

			double Q = (Qub + Qlb) / 2;
			do {
				double gradient = this.simulateCycle(Q, inventory, reducedHorizon)[0]
						- this.simulateCycle(Q + step, inventory, reducedHorizon)[0];
				if (gradient > 0) {
					Qlb = Q + step;
				} else {
					Qub = Q;
				}
				Q = (Qub + Qlb) / 2;
			} while ((reducedHorizon[t] instanceof DiscreteDistributionInt && Qub - Qlb > step));
			Q = Math.round(Q);

			double curCapital = (this.simulateCycle(Q, inventory, reducedHorizon)[0]) / reducedHorizon.length;
			if (curCapital < bestCapital) {
				bestCapital = curCapital;
				bestQ = Q;
			}
		}

		return noOrder < bestCapital ? 0 : bestQ;
	}

	public double simulateSingleRun() {
		double i = this.inventory;
		double[][] d = Sampling.generateLHSamples(distributions, 1);
		double c = 0;
		for (int t = 0; t < this.distributions.length; t++) {
			double Q = getQ(t, i);
			if (Q > 0) {
				c += orderingCost(Q);
				i += Q;
			}

			i -= d[t][0];
			c += inventoryCost(i);
		}
		return c;
	}

	public double[] simulate() {
		Sampling.resetNextSubstream();
		TallyStore observationsTally = new TallyStore();
		for (int runs = 0; runs < this.minMonteCarloSimulationRunsP; runs++) {
			observationsTally.add(simulateSingleRun());
		}
		double[] centerAndRadius = new double[2];
		observationsTally.confidenceIntervalStudent(confidenceLevel, centerAndRadius);
		/*
		 * while(Math.abs(centerAndRadius[1]/centerAndRadius[0]) >=
		 * percentageError/100){ observationsTally.add(simulateSingleRun());
		 * observationsTally.confidenceIntervalStudent(confidenceLevel,
		 * centerAndRadius); }
		 */
		/*
		 * IntStream.iterate(0, i -> i + 1) .limit(observationsTally.numberObs())
		 * .forEach(i -> System.out.println(observationsTally.getArray()[i]));
		 */
		return centerAndRadius;
	}

	public static void runSingleInstance() {
		double K = 500;
		double v = 2;
		double pai = 20;
		double h = 1;
		double maxOrderQuantity = 100;

		double lambda[] = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20 };
		Distribution[] distributions = IntStream.iterate(0, i -> i + 1).limit(lambda.length)
				.mapToObj(i -> new PoissonDist(lambda[i])).toArray(Distribution[]::new);

		double inventory = 0;
		RRSimOpt a = new RRSimOpt(K, v, pai, h, distributions, inventory, maxOrderQuantity);

		double[] stats = a.simulate();
		System.out.println(stats[0] + " " + stats[1]);
	}

	public static void main(String args[]) {
		long currTime1 = System.currentTimeMillis();
		runSingleInstance();
		double time = (System.currentTimeMillis() - currTime1) / 1000;
		System.out.println("running time is " + time + " s");
	}
}