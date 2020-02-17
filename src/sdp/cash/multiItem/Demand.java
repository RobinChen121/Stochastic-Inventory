/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jan 20, 2020, 2:14:49 PM
 * @Desc: 
 *
 *
 * 
 */
package sdp.cash.multiItem;

/**
 * a class for random demand of one product in two items, 
 * demand are assumed to be integers
 */
public class Demand {
	int demand;
	
	public Demand(int demand) {
		this.demand = demand;
	}
	
	public int getDemand() {
		return demand;
	}	

}