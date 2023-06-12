package sdp.cash;

import sdp.inventory.State;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: June 10, 2023---3:22:49 PM
*@description:  states for the cash survival problem, including 3 states
*/

public class RiskState extends CashState{
	boolean bankruptBefore = false;

	public RiskState(int period, double initialInventory, double iniCash, boolean bankruptBefore) {
		super(period, initialInventory, iniCash);
		this.bankruptBefore = false;
	}
	
	public double getIniCash() {
		return this.iniCash;
	}
	
	public boolean getBankruptBefore() {
		return this.bankruptBefore;
	}
	
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + initialInventory + iniCash + bankruptBefore;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof RiskState)
			return ((RiskState) o).period == this.period &&
					((RiskState) o).initialInventory == this.initialInventory &&
							((RiskState) o).iniCash == this.iniCash &&
									((RiskState) o).bankruptBefore == this.bankruptBefore;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory = " + initialInventory + ", iniCash = " + iniCash
								+ ", bankuptBefore = " + bankruptBefore; 
	}
	
}

