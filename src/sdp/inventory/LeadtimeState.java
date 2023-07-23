package sdp.inventory;
import sdp.inventory.State;

/**
*@author: zhenchen
*@date: Jul 23, 2023, 10:34:26 AM
*@desp: TODO
*
*/

public class LeadtimeState extends State{
	double preInventory;
	
	public LeadtimeState(int period, double initialInventory, double preInventory) {
		super(period, initialInventory);
		this.preInventory = preInventory;
	}

}


