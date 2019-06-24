/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Jun 21, 2019, 1:20:06 PM
 * @Desc: 
 *
 *
 * 
 */
package multiItem;

public interface ImmediateValueFunction<S, A1, A2, R, D> {
	
	/**
	 * The immediate value function
	 * 
	 * @param initialState the initial state of the stochastic process.
	 * @param action1 the chosen action1.
	 * @param action2 the chosen action2.
	 * @param randomOutcome the final state of the stochastic process.
	 * @return the immediate value of a transition from {@code initialState} to {@code finalState} under twp chosen
	 * actions  {@code action1}, {@code action2}.
	 */
	public D apply (S initialState, A1 action1, A2 action2, R randomOutcome);

}
