package sdp.inventory;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 12, 2018---11:47:12 AM
*@description:  a state transition interface
*/

public class StateTransition {
	
	/**
	 * 
	 * @param <S> initial state
	 * @param <A> action
	 * @param <R> random demand
	 * @param <S2> output state
	 */
	
	public interface StateTransitionFunction<S, A, R, S2> {
		public S2 apply(S s, A a, R r);
	}
	
	
	public interface StateTransitionFunctionV <S, R, S2>{
		public S2 apply (S s, R r);
	}
	

}
