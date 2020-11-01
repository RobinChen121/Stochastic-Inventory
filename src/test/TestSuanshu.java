/**
 * @date: Nov 1, 2020
 */
package test;

import java.awt.print.Printable;

import com.numericalmethod.suanshu.vector.doubles.Vector;
import com.numericalmethod.suanshu.vector.doubles.dense.DenseVector;

/**
 * @author: Zhen Chen
 * @email: 15011074486@163.com
 * @date: Nov 1, 2020
 * @Desc:
 *
 */
public class TestSuanshu {


	public static void main(String[] args) {
		Vector v1 = new DenseVector(new double[] {1, 2, 3, 4, 5});
	    Vector v2 = new DenseVector(new double[] {5, 4, 3, 2, 1});
	    Vector v3 = v1.add(v2);
	    //log.info("Adding vectors: {}", v3);
	    System.out.println(v3);
	}

}
