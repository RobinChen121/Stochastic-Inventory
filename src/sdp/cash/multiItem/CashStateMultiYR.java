/**
 * @date: Jul 7, 2020
 */
package sdp.cash.multiItem;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 7, 2020
 * @Desc:
 *
 */
public class CashStateMultiYR {
	int period;
	double iniInventory1;
	double iniInventory2;
	double iniR;
	
	public CashStateMultiYR(int t, double iniInventory1, double iniInventory2, double R) {
		this.period = t;
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
		if (o instanceof CashStateMultiYR)
			return ((CashStateMultiYR) o).period == this.period &&
					((CashStateMultiYR) o).iniInventory1 == this.iniInventory1 &&
							((CashStateMultiYR) o).iniInventory2 == this.iniInventory2 &&
								(int) ((CashStateMultiYR) o).iniR == (int) this.iniR;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory1 = " + iniInventory1 + ", iniInventory2 = " + iniInventory2 +  ", iniR = " + iniR; 
	}

}
