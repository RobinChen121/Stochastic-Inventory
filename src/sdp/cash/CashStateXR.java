package sdp.cash;

import sdp.inventory.State;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Feb 13, 2020---12:35:49 PM
*@description:  State class for cash flow problem in Chao (2008), states are initial inventory x and working capital 
*               R = w + cx
*/

public class CashStateXR extends State{
	double iniR = 0;
	double unitVariCost = 0;

	public CashStateXR(int period, double iniInventory, double R, double variCost) {
		super(period, iniInventory);
		this.iniR = R;
		this.unitVariCost = variCost;
	}
	
	public double getIniInventory() {
		return this.initialInventory;
	}
	
	
	public double getIniR() {
		return this.iniR;
	}
	
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + initialInventory + iniR;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof CashStateXR)
			return ((CashStateXR) o).period == this.period &&
					((CashStateXR) o).initialInventory == this.initialInventory &&
							((CashStateXR) o).iniR == this.iniR;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory = " + initialInventory + ", iniR = " + iniR; 
	}
	
}
