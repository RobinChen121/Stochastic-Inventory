/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 21, 2019, 11:16:49 AM
 * @Desc:  cash state for multi item: x1, x2, w
 *
 *
 * 
 */
package sdp.cash.multiItem;



public class CashStateMulti {
	int period;
	int iniInventory1;
	int iniInventory2;
	double iniCash;
	
	public CashStateMulti(int period, int iniInventory1, int iniInventory2, double iniCash) {
		this.period = period;
		this.iniInventory1 = iniInventory1;
		this.iniInventory2 = iniInventory2;
		this.iniCash = iniCash;		
	}
	
	public int getPeriod() {
		return this.period;
	}
	
	public int getIniInventory1() {
		return this.iniInventory1;
	}
	
	public int getIniInventory2() {
		return this.iniInventory2;
	}
	
	public double getIniCash() {
		return this.iniCash;
	}
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + iniInventory1 + iniInventory2 + iniCash + iniCash;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof CashStateMulti)
			return ((CashStateMulti) o).period == this.period &&
					((CashStateMulti) o).iniInventory1 == this.iniInventory1 &&
							((CashStateMulti) o).iniInventory2 == this.iniInventory2 &&
								(int) ((CashStateMulti) o).iniCash == (int) this.iniCash;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory1 = " + iniInventory1 + ", iniInventory2 = " + iniInventory2 +  ", iniCash = " + iniCash; 
	}

}
