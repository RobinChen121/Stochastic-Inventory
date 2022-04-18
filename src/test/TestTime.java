/**
 * @date: Apr 23, 2020
 */
package test;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Apr 23, 2020
 * @Desc: output running time in float seconds
 *
 */
public class TestTime {

	public static void main(String[] args) {
		long currTime = System.currentTimeMillis();
		double sum = 0;
		for (int i  = 0; i < 10000000; i++)
			sum += i;
		
		double time = (System.currentTimeMillis() - currTime) / 1000.0; //the denominator should be float in order to output float number
		System.out.println("running time is " + time + "s");

	}

}
