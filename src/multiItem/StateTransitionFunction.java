/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 24, 2019, 10:09:26 AM
 * @Desc: 
 *
 *
 * 
 */
package multiItem;


	
@FunctionalInterface
public interface StateTransitionFunction <S, A, R, D> { 
	/**
	 * The state transition function
	 * 
	 * @param initialState the initial state of the stochastic process.
	 * @param actions the chosen actions.
	 * @param RandomDemands random demands.
	 * @return the final state given {@code initialState} under chosen {@code actions} and {@code randomDemands}.
	 */
	public D apply (S initialState, A actions, R randomDemands);
}

