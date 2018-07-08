package forTest;

public class TestClass {
	int k;
	int g;
	
	public TestClass(int k) {
		this.k = k;
		this.g = getG(k);
	}
	
	void outPut() {
		System.out.println(k);
		System.out.println(g);
	}
	
	int getG(int k) {
		return k+1;
	}
	
	public static void main(String[] args) {
		TestClass temp = new TestClass(5);
		temp.outPut();
		
	}

}
