package sdp.inventory;


public class CheckKConvexity {
		
	public static String checkCK(double[][] yG, double fixOrderCost, int capacity) {
		double minInventorys = yG[0][0];
		double maxInventorys = yG[yG.length - 1][0]; 
		int xLength = (int) (maxInventorys - minInventorys + 1);
		
		int mark = 0;
		// follow the definition of CK convex in Gallego and Scheller-Wolf (2000)
		for (int y = 0; y < xLength; y++)  // y
			for (int z = 0; z < capacity; z++)
				for (int b = 1; b < capacity ; b++) { // c is at least 1 for GB(y)
					if (y - b <= 0 || y + z >= xLength || yG[y + z][1] + fixOrderCost > yG[y][1] + z * (yG[y][1] - yG[y - b][1])/b - 0.1)
						continue;
					else {
						System.out.println(yG[y + z][1] + fixOrderCost);
						System.out.println(yG[y][1] + z * (yG[y][1] - yG[y - b][1])/b);
						int xy = (int) yG[y][0];
						int xb = (int) yG[b][0];
						int xz = (int) yG[z][0];
						System.out.printf("z = %d, y = %d, b = %d", xz, xy, xb);
						System.out.println();
						System.out.println("not CK convex");
						return "not CK convex";
					}
				}
		
		if (mark == 0) {
			System.out.println("CK convexity holds");
			return "CK convexity holds";
		}
		return null;
	}
	
	
	public static String check(double[][] yG, double fixOrderCost) {
		double minInventorys = yG[0][0];
		double maxInventorys = yG[yG.length - 1][0]; 
		int xLength = (int) (maxInventorys - minInventorys + 1);
		
		int mark = 0;
		for (int a = 0; a < xLength; a++)
			for (int b = 0; b < a; b++)
				for (int c = 0; c < b ; c++) {
					if (b - c == 0 || yG[a][1] + fixOrderCost > yG[b][1] + (a - b) * (yG[b][1] - yG[c][1])/(b - c) - 0.1)
						continue;
					else {
						System.out.println(yG[a][1] + fixOrderCost);
						System.out.println(yG[b][1] + (a - b) * (yG[b][1] - yG[c][1])/(b - c));
						int xb = (int) yG[c][0];
						int x = (int) yG[b][0];
						int xa = (int) yG[a][0];
						System.out.printf("x-b = %d, x = %d, x+a = %d", xb, x, xa);
						System.out.println();
						System.out.println("not K convex");
						return "not K convex";
					}
				}
		
		if (mark == 0) {
			System.out.println("CK convexity holds");
			return "CK convexity holds";
		}
		return null;
	}
}
