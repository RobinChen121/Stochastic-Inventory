package sdp.inventory;

/**
*@author: zhenchen
*@date: Jul 23, 2023, 10:34:26 AM
*@desp: TODO
*
*/

public class LeadtimeState extends State{
	double preQ;
	
	public LeadtimeState(int period, double initialInventory, double preQ) {
		super(period, initialInventory);
		this.preQ = preQ;
	}
	
	public double getPreQ() {
		return this.preQ;
	}
	
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + initialInventory + preQ;
		return hash.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof LeadtimeState)
			return ((LeadtimeState) o).period == this.period &&
					((LeadtimeState) o).initialInventory == this.initialInventory &&
							((LeadtimeState) o).preQ == this.preQ;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period + ", " + "initialInventory = " + initialInventory + ", " + "preQ = " + preQ;
	}
	
	
	public int compareTo(LeadtimeState o) {
		return this.period > o.period ? 1
				: this.period == o.period
						? this.initialInventory > o.initialInventory ? 1
								: this.initialInventory == o.initialInventory ? 
										this.preQ > o.preQ ?  1:
											(this.preQ == o.preQ ? 0 : -1):
												-1:-1;
	}
		

}


