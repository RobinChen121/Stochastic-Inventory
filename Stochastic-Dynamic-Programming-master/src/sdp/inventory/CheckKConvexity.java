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
					if (yG[a][1] + fixOrderCost > yG[b][1] + (a - b)/(b - c) * (yG[b][1] - yG[c][1]) - 0.1)
						continue;
					else {
						System.out.print("not K convext");
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
