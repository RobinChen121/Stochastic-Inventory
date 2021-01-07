package sdp.inventory;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 12, 2018---12:06:39 PM
*@description:  immediate value function for stochastic dynamic programming
*/

public class ImmediateValue {
	
	/**
	 * 
	 * @author Administrator
	 *
	 * @param <S> state
	 * @param <A> action
	 * @param <R> random demand
	 * @param <V> output value
	 */
	

	public interface ImmediateValueFunction<S, A, R, V> {
		public V apply(S s, A a, R r);
	} 
	
	

	public interface ImmediateValueFunctionV <S, R, V>{
		public V apply (S s, R r);
	}


}
