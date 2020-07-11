/**
 * @date: Jul 10, 2020
 */
package sdp.cash.multiItem;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jul 10, 2020
 * @Desc: for computing y*(R)
 *
 */
public class CashStateR {
	int period;

	double iniR;
	
	public CashStateR(int t,double R) {
		this.period = t;
		this.iniR = R;		
	}
	
	public int getPeriod() {
		return this.period;
	}
	
	
	public double getIniR() {
		return this.iniR;
	}
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + iniR;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof CashStateMultiYR)
			return ((CashStateMultiYR) o).period == this.period &&
								(int) ((CashStateMultiYR) o).iniR == (int) this.iniR;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period + ", iniR = " + iniR; 
	}

}
