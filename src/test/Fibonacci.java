package test;

/**
*@author: zhen chen
*@date: Feb 22, 2025, 10:00:44 PM
*@desp: TODO
*
*/

import java.util.HashMap;
import java.util.Map;
import java.math.BigInteger;

public class Fibonacci {
	private static Map<Integer, BigInteger> cache = new HashMap<>();

	// 递归实现（带缓存优化）
	public static BigInteger fibonacciMemoization(int n) {
		if (n <= 1) {
			return BigInteger.valueOf(n);
		}

		// 检查缓存中是否已经计算过
		if (cache.containsKey(n)) {
			return cache.get(n);
		}

		// 计算并缓存结果
		BigInteger result = fibonacciMemoization(n - 1).add(fibonacciMemoization(n - 2));
		cache.put(n, result);

		return result;
	}

	public static void main(String[] args) {
		int n = 2500;
		long currTime = System.currentTimeMillis();
		System.out.println("Fibonacci(" + n + ") = " + fibonacciMemoization(n));
		double time = (System.currentTimeMillis() - currTime)/ 1000.0;
		String formatted = String.format("%.4f", time); 
		System.out.println("running time is " + formatted + "s");
	}
}
