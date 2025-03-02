package test;

/**
 * @author: zhenchen
 * @date: Feb 22, 2025, 10:28:36 PM
 * @desp: TODO
 */

public class Fibonacci2 {
    // Recursive Fibonacci
    public static long fibRecursive(int n) {
        if (n <= 1) {
            return n;
        } else {
            return fibRecursive(n - 1) + fibRecursive(n - 2);
        }
    }

    // Iterative Fibonacci
    public static long fibIterative(int n) {
        long a = 0, b = 1;
        for (int i = 0; i < n; i++) {
            long temp = a + b;
            a = b;
            b = temp;
        }
        return a;
    }

    public static void main(String[] args) {
        long startTime = System.nanoTime();
        fibRecursive(40);
        long endTime = System.nanoTime();
        System.out.println("Java recursive Fibonacci time: " + (endTime - startTime) / 1000000000.0 + " seconds");

        startTime = System.nanoTime();
        fibIterative(30);
        endTime = System.nanoTime();
        System.out.println("Java iterative Fibonacci time: " + (endTime - startTime) / 1000000000.0 + " seconds");
    }
}


