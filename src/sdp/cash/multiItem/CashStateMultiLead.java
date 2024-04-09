package sdp.cash.multiItem;

/**
*@author: zhenchen
*@date: Apr 2, 2024, 7:56:34 PM
*@desp: TODO
*
*/

public class CashStateMultiLead {
	int period;
	double iniInventory1;
	double iniInventory2;
	double iniCash;
	double preQ1;
	double preQ2;
	
	public CashStateMultiLead(int period, double iniInventory1, double iniInventory2, double preQ1, double preQ2, double iniCash) {
		this.period = period;
		this.iniInventory1 = iniInventory1;
		this.iniInventory2 = iniInventory2;
		this.iniCash = iniCash;		
		this.preQ1 = preQ1;
		this.preQ2 = preQ2;
	}
	
	public int getPeriod() {
		return this.period;
	}
	
	public double getIniInventory1() {
		return this.iniInventory1;
	}
	
	public double getIniInventory2() {
		return this.iniInventory2;
	}
	
	public double getPreQ1() {
		return this.preQ1;
	}
	
	public double getPreQ2() {
		return this.preQ2;
	}
	
	public double getIniCash() {
		return this.iniCash;
	}
	
	public double getIniR(double[] variCosts) {
		return this.iniCash + variCosts[0] * this.iniInventory1 + variCosts[1] * this.iniInventory2;
	}
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + iniInventory1 + iniInventory2 + iniCash + iniCash;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof CashStateMultiLead)
			return ((CashStateMultiLead) o).period == this.period &&
					((CashStateMultiLead) o).iniInventory1 == this.iniInventory1 &&
							((CashStateMultiLead) o).iniInventory2 == this.iniInventory2 &&
									((CashStateMultiLead) o).preQ1 == this.preQ1 &&
									((CashStateMultiLead) o).iniInventory2 == this.iniInventory2 &&
										(int) ((CashStateMultiLead) o).iniCash == (int) this.iniCash;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory1 = " + iniInventory1 + ", iniInventory2 = " + iniInventory2 +", "+"preQ1 = " + preQ1 + ", preQ2 = " + preQ2+  ", iniCash = " + iniCash; 
	}

}


