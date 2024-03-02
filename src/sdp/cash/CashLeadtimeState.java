package sdp.cash;


/**
*@author: zhenchen
*@date: Mar 1, 2024, 7:01:49 PM
*@desp: TODO
*
*/

public class CashLeadtimeState extends CashState{
	double preQ;
	
	public CashLeadtimeState(int period, double initialInventory, double iniCash, double preQ) {
		super(period, initialInventory, iniCash);
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
		if (o instanceof CashLeadtimeState)
			return ((CashLeadtimeState) o).period == this.period &&
					((CashLeadtimeState) o).initialInventory == this.initialInventory &&
							((CashLeadtimeState) o).iniCash == this.iniCash &&
								((CashLeadtimeState) o).preQ == this.preQ;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period + ", " + "initialInventory = " + initialInventory + ", " + "iniCash = " + iniCash + ", "+ "preQ = " + preQ;
	}

}


