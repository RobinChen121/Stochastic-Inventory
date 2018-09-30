package sdp.inventory;


public class CheckKConvexity {
	
	public static String check(double[][] yG, double fixOrderCost) {
		double minInventorys = yG[0][0];
		double maxInventorys = yG[yG.length - 1][0]; 
		int xLength = (int) (maxInventorys - minInventorys + 1);
		
		int mark = 0;
		for (int a = 0; a < xLength; a++)
			for (int b = 0; b < a; b++)
				for (int c = 0; c < b; c++) {
					if (yG[a][1] + fixOrderCost > yG[b][1] + (a - b) * (yG[b][1] - yG[c][1])/(b - c) - 5)
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
			System.out.println("K convexity holds");
			return "K convexity holds";
		}
		return null;
	}

}
