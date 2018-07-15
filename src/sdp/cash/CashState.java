package sdp.cash;

import sdp.inventory.State;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 13, 2018---3:35:49 PM
*@description:  State calss for cash flow problem
*/

public class CashState extends State{
	double iniCash;

	public CashState(int period, double initialInventory, double iniCash) {
		super(period, initialInventory);
		this.iniCash = iniCash;
	}
	
	public double getIniCash() {
		return this.iniCash;
	}
	
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + initialInventory + iniCash;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof CashState)
			return ((CashState) o).period == this.period &&
					((CashState) o).initialInventory == this.initialInventory &&
							((CashState) o).iniCash == this.iniCash;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory = " + initialInventory + ", iniCash = " + iniCash; 
	}
	
}
