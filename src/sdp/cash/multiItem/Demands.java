/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 24, 2019, 11:20:24 AM
 * @Desc: 
 *
 *
 * 
 */
package sdp.cash.multiItem;



/**
 * a class for random demand of two items, 
 * demand are assumed to be integers
 */
public class Demands {
	int demand1;
	int demand2;
	
	public Demands(int demand1, int demand2) {
		this.demand1 = demand1;
		this.demand2 = demand2;
	}
	
	public int getFirstDemand() {
		return demand1;
	}
	
	public int getSecondDemand() {
		return demand2;
	}

}
