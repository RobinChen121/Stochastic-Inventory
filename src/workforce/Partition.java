package workforce;

public class Partition {

    public static void partition1(int n) {
        partition1(n, n, "");
    }
    public static void partition1(int n, int max, String prefix) {
        if (n == 0) {
            System.out.println(prefix);
            return;
        }

        for (int i = Math.min(max, n); i >= 1; i--) {
            partition1(n-i, i, prefix + " " + i);
        }
    }
    
	public static void partion(int n, int k) {
		partion(n, k, "");	
	}
	
	public static void partion(int n, int k, String str) {		
		if (k == 0)
			return;
		if (k == 1) {
			System.out.println(str + " " + n);
			return;
		}
		
		for (int i = 0; i <= n; i++) {
			partion(n-i, k-1, str + " " + i);
		}
	}

    public static void main(String[] args) {
        int n = 6;
        int k = 3;
        partion(n, k);
        //partition1(n);
    }

}