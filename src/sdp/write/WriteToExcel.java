package sdp.write;

import java.io.FileWriter;
import java.io.IOException;


/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 9, 2018---12:48:56 PM
*@Description:  read an array to txt or excel files
*/

public class WriteToExcel {
	
	public void writeArrayToTxt(double[][] data, String string) {
		int rowNum = data.length;
		int columnNum = data[0].length;
		try {
			FileWriter fw = new FileWriter(string);
			for (int i = 0; i < rowNum; i++) {
				for (int j = 0; j < columnNum; j++)
					fw.write(data[i][j]+ "\t");
				fw.write("\n");
			}
			fw.close();
		}
		catch (IOException e){
			e.printStackTrace();
		}		
	}
	
	// 导出到 excel的代码其实跟导出到 txt 的代码一样
	public void writeArrayToExcel(double[][] data, String string) {
		int rowNum = data.length;
		int columnNum = data[0].length;
		try {
			FileWriter fw = new FileWriter(string);
			for (int i = 0; i < rowNum; i++) {
				for (int j = 0; j < columnNum; j++)
					fw.write(data[i][j]+ "\t");
				fw.write("\n");
			}
			fw.close();
		}
		catch (IOException e){
			e.printStackTrace();
		}		
	}
	
	
	/** output one-dimension array in excel
	 * @param data
	 * @param string
	 * @date: May 31, 2020, 6:05:18 PM 
	 */
	public void writeToExcelAppend(double[] data, String string) {
		int columnNum = data.length;
		try {
			FileWriter fw = new FileWriter(string, true);
			for (int i = 0; i < columnNum; i++) {
				fw.write(data[i]+ "\t");			
			}
			fw.write("\n");
			fw.close();
		}
		catch (IOException e){
			e.printStackTrace();
		}		
	}
	
	
	public void writeArrayToExcel(double[][] data, String string, String head) {
		int rowNum = data.length + 1;
		int columnNum = data[0].length;
		try {
			FileWriter fw = new FileWriter(string);
			for (int i = 0; i < rowNum; i++) {
				if (i == 0) {
					fw.write(head);
					fw.write("\n");
				}
				else {
					for (int j = 0; j < columnNum; j++)
						fw.write(data[i - 1][j] + "\t");
					fw.write("\n");
				}
			}
			fw.close();
		}
		catch (IOException e){
			e.printStackTrace();
		}
		
	}
	
	
	public static void main(String[] args) {

		 double[][] demands = {{7,7,7,7,7,7},
                 {2,3,4,5,6,7},
                 {8,7,6,5,4,3},
                 {5,6,7,8,7,6},
                 {8,5,2,1,2,5},
                 {8,4,1,3,1,3},
                 {1,3,8,4,8,7},
                 {1,4,7,3,5,8},
                 {3,8,4,4,6,2},
                 {3,1,5,8,4,4}
                 };
		 
		 WriteToExcel wa = new WriteToExcel();
		 String headStrings = "c1" + "\t" + "c2" + "\t" + "c3" + "\t" + "c4" + "\t" + "c5" + "\t" + "c6";
		 //wa.writeArrayToTxt(demands, "mytxt.txt");
		 wa.writeArrayToExcel(demands, "mytxt.xls", headStrings);
	}
}
