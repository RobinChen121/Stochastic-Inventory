/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 24, 2019, 10:34:48 AM
 * @Desc: 
 *
 *
 * 
 */
package sdp.cash.multiItem;



/**
 * a class for actions in the double item SDP
 */
public class Actions {
	int firstAction;
	int secondAction;
	
	public Actions(int action1, int action2) {
		this.firstAction = action1;
		this.secondAction = action2;		
	}
	
	public int getFirstAction() {
		return firstAction;
	}
	
	public int getSecondAction() {
		return secondAction;
	}

}
