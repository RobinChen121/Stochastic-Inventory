package sdp.inventory;


import sdp.inventory.Recursion.OptDirection;

/**
* @author Zhen Chen
* @date: 2018年12月12日 下午4:14:52  
* @email: 15011074486@163.com,
* @licence: MIT licence. 
*
* @Description:  record best action and its corresponding expected value for a feasible action
*/

public class BestActionValue {
	double bestAction = 0;
	double bestValue;
	OptDirection direction;
	
	public BestActionValue(OptDirection direction) {
		this.direction = direction;
		bestValue = direction == OptDirection.MAX ? -Double.MAX_VALUE : Double.MAX_VALUE;
	}
	
	public double getBestAction() {
		return bestAction;
	}
	
	public double getBestValue() {
		return bestValue;
	}
	
	public synchronized void update(double action, double value) {
		if (direction == OptDirection.MIN) {
			if (value < bestValue) {
				bestValue = value;
				bestAction = action;
			}				
		}
		else {
			if (value > bestValue) {
				bestValue = value;
				bestAction = action;
			}
		}
			
	}

}


