package sdp.cash;

import sdp.inventory.State;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 13, 2018---3:35:49 PM
*@description:  State class for cash flow problem, cash can be integers(seems not influence much to the computation time)
*/

public class CashState extends State{
	double iniCash = 0;
	int iniCashInteger = 0;

	public CashState(int period, double initialInventory, double iniCash) {
		super(period, initialInventory);
		this.iniCash = iniCash;
	}
	
	public CashState(int period, double initialInventory, int iniCash) {
		super(period, initialInventory);
		this.iniCash = iniCash;
	}
	
	public double getIniCash() {
		if (iniCash > 0)
			return this.iniCash;
		else {
			return this.iniCashInteger;
		}
	}
	
	public double getIniCashInteger() {
		return this.iniCashInteger;
	}
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + initialInventory + iniCash + iniCashInteger;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof CashState)
			return ((CashState) o).period == this.period &&
					((CashState) o).initialInventory == this.initialInventory &&
							((CashState) o).iniCash == this.iniCash &&
									((CashState) o).iniCashInteger == this.iniCashInteger;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory = " + initialInventory + ", iniCash = " + iniCash + ", iniCashInteger = " + iniCashInteger; 
	}
	
}
