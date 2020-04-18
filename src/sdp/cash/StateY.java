/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Apr 17, 2020, 6:43:09 PM
 * @Desc: state for G()
 *
 *
 * 
 */
package sdp.cash;


public class StateY {
	protected int period;
	protected double iniY;
	
	public StateY(int period, double y) {
		this.period = period;		
		this.iniY = y;		
	}
	
	public int getPeriod() {
		return this.period;
	}

	public double getIniY() {
		return this.iniY;
	}
	
	@Override
	public int hashCode() {
		String hash = "";
		hash = hash + period + iniY;
		return hash.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof StateY)
			return ((StateY) o).period == this.period && ((StateY) o).iniY == this.iniY;
		else
			return false;
	}

	@Override
	public String toString() {
		return "period = " + period + ", " + "initialInventory = " + iniY;
	}

	public int compareTo(StateY o) {
		return this.period > o.period ? 1
				: this.period == o.period
						? this.iniY > o.iniY ? 1
								: (this.iniY == o.iniY ? 0 : -1)
						: -1;
	}

}
