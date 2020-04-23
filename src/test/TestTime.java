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
		
		double time = (System.currentTimeMillis() - currTime) / 1000.0; // 分母改为小数就输出的为小数，否则输出为整数
		System.out.println("running time is " + time + "s");

	}

}
