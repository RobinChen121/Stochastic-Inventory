package sdp.cash.multiItem;

/**
 * @author chen
 * @email: 15011074486@163.com
 * @Date: 2021 Feb 26, 21:36:35
 * @Description: TODO 
 * 
 */
public class CashStateMultiXYW {
	int period;
	double iniInventory1;
	double iniInventory2;
	double y1;
	double y2;
	double iniW;
	
	public CashStateMultiXYW(int t, double iniInventory1, double iniInventory2, double y1, double y2, double w) {
		this.period = t;
		this.iniInventory1 = iniInventory1;
		this.iniInventory2 = iniInventory2;
		this.y1 = y1;
		this.y2 = y2;
		this.iniW = w;		
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
	
	public double getY1() {
		return this.y1;
	}
	
	public double getY2() {
		return this.y2;
	}
	
	public double getIniW() {
		return this.iniW;
	}
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + iniInventory1 + iniInventory2 + y1 + y2 + iniW;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof CashStateMultiXYW)
			return ((CashStateMultiXYW) o).period == this.period &&
					((CashStateMultiXYW) o).iniInventory1 == this.iniInventory1 &&
							((CashStateMultiXYW) o).iniInventory2 == this.iniInventory2 &&
									((CashStateMultiXYW) o).y1 == this.y1 &&
										((CashStateMultiXYW) o).y2 == this.y2 &&
											(int) ((CashStateMultiXYW) o).iniW == (int) this.iniW;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"iniInventory1 = " + iniInventory1 + ", iniInventory2 = " + iniInventory2 + ", y1= "+y1 +", y2 = " + y2+", iniW = " + iniW; 
	}

}
