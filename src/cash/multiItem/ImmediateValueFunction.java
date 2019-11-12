/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 21, 2019, 1:20:06 PM
 * @Desc: 
 *
 *
 * 
 */
package cash.multiItem;

/**
 * The immediate value function
 * 
 * @param initialState the initial state of the stochastic process.
 * @param actions the chosen actions.
 * @param randomDemands the random demand of the stochastic process.
 * @return the immediate value of a transition from {@code initialState} under two chosen
 * actions  {@code actions} and stochastic demands {@code randomDemands}.
 */

public interface ImmediateValueFunction<S, A, R, D> {
	
	
	public D apply (S initialState, A actions, R randomDemands);

}
