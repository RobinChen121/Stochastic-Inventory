/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 21, 2019, 11:16:49 AM
 * @Desc:  cash state for multi item: x1, x2, R
 *
 *
 * 
 */
package sdp.cash.multiItem;



public class CashStateMultiXR {
	int period;
	double iniInventory1;
	double iniInventory2;
	double iniR;
	
	public CashStateMultiXR(int d, double iniInventory1, double iniInventory2, double R) {
		this.period = d;
		this.iniInventory1 = iniInventory1;
		this.iniInventory2 = iniInventory2;
		this.iniR = R;		
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
	
	public double getIniR() {
		return this.iniR;
	}
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + iniInventory1 + iniInventory2 + iniR;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof CashStateMultiXR)
			return ((CashStateMultiXR) o).period == this.period &&
					((CashStateMultiXR) o).iniInventory1 == this.iniInventory1 &&
							((CashStateMultiXR) o).iniInventory2 == this.iniInventory2 &&
								(int) ((CashStateMultiXR) o).iniR == (int) this.iniR;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory1 = " + iniInventory1 + ", iniInventory2 = " + iniInventory2 +  ", iniR = " + iniR; 
	}

}
