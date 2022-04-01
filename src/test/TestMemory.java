/**
 * @date: May 17, 2020
 */
package test;

import java.util.stream.IntStream;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: May 17, 2020
 * @Desc: test java memory limit
 *
 */
public class TestMemory {

	/**
	 * @param args
	 * @date: May 17, 2020, 11:37:49 PM 
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
//		int[][] test = new int[122281250][3];
//		System.out.println("this test");
		
		int[] a = {1, 2, 3, 4};
		int b = IntStream.of(a).skip(4).reduce(1, (x, y) -> x*y);
		System.out.println(b);

	}

}

