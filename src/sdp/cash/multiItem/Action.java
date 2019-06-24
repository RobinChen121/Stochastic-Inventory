/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 19, 2019, 9:33:24 PM
 * @Desc: 
 *
 *
 * 
 */
package sdp.cash.multiItem;

import java.io.Serializable;

public abstract class Action implements Serializable {
	
	private static final long serialVersionUID = 1L;
	
	public abstract boolean equals(Object Action);
	public abstract int hashcode();

}
