/**
 * @date: Apr 28, 2020
 */
package test;

import java.util.ArrayList;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Apr 28, 2020
 * @Desc:
 *
 */
public class TestArrayListDoubleArray {

	/**
	 * @param args
	 * @date: Apr 28, 2020, 10:31:23 PM 
	 */
	public static void main(String[] args) {
		ArrayList<double[]> list = new ArrayList<>();
		list.add(new double[] {1, 2});
		System.out.println(list.get(0)[0]);

	}

}
