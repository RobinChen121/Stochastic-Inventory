/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Apr 17, 2020, 6:43:09 PM
 * @Desc: state for period to record the bestY in each period
 *
 *
 * 
 */
package sdp.cash;


public class StateP {
	protected int period;
	
	public StateP(int period) {
		this.period = period;			
	}
	
	public int getPeriod() {
		return this.period;
	}


	
	@Override
	public int hashCode() {
		String hash = "";
		hash = hash + period;
		return hash.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof StateP)
			return ((StateP) o).period == this.period;
		else
			return false;
	}

	@Override
	public String toString() {
		return "period = " + period ;
	}

	public int compareTo(StateP o) {
		return this.period > o.period ? 1
				: this.period == o.period
						? 0 : -1;
	}

}
