package sdp.inventory;


/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 9, 2018---10:24:14 PM
 * @description: A State class with inventory as state for dynamic programming
 */

public class State implements Comparable<State> {
	protected int period;
	protected double initialInventory;

	public State(int period, double initialInventory) {
		this.period = period;
		this.initialInventory = initialInventory;
	}

	
	public int getPeriod() {
		return this.period;
	}

	public double getIniInventory() {
		return this.initialInventory;
	}

//	public double[] getFeasibleActions(double maxOrderQuantity, double stepSize) {
//		double[] feasibleActions = new double[(int) (maxOrderQuantity/stepSize) + 1];
//		int index = 0;
//		for (double i = 0; i <= maxOrderQuantity; i = i + stepSize) {
//			feasibleActions[index] = i;
//			index++;
//		}
//		return feasibleActions;
//	}
//
//	public State stateTransition(double action, double randomDemand, double minInventory,
//			double maxInventory, boolean isForDrawGy) {
//		double nextInventory = isForDrawGy && this.period == 1 ? this.initialInventory - randomDemand
//				: this.initialInventory + action - randomDemand;
//		nextInventory = nextInventory > maxInventory ? maxInventory : nextInventory;
//		nextInventory = nextInventory < minInventory ? minInventory : nextInventory;
//		return new State(this.period + 1, nextInventory);
//	}
//	
//	public State stateTransition(double action, double randomDemand, double minInventory,
//			double maxInventory) {
//		double nextInventory = this.initialInventory + action - randomDemand;
//		nextInventory = nextInventory > maxInventory ? maxInventory : nextInventory;
//		nextInventory = nextInventory < minInventory ? minInventory : nextInventory;
//		return new State(this.period + 1, nextInventory);
//	}
	

	@Override
	public int hashCode() {
		String hash = "";
		hash = hash + period + initialInventory;
		return hash.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof State)
			return ((State) o).period == this.period && ((State) o).initialInventory == this.initialInventory;
		else
			return false;
	}

	@Override
	public String toString() {
		return "period = " + period + ", " + "initialInventory = " + initialInventory;
	}

	@Override
	public int compareTo(State o) {
		return this.period > o.period ? 1
				: this.period == o.period
						? this.initialInventory > o.initialInventory ? 1
								: (this.initialInventory == o.initialInventory ? 0 : -1)
						: -1;
	}
}
