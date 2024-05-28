package test;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;

/**
*@author: zhenchen
*@date: May 25, 2024, 10:24:58 AM
*@desp: TODO
*
*/

public class TestDist {
	public static void main(String[] args) {
		Distribution dist = new NormalDist();
		System.out.println(dist.getMean());
	}

}


