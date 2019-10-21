package sdp.write;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;

/**
*@author: Zhen Chen
*@email: 15011074486@163.com
*@date: Jul 12, 2018---10:00:09 PM
*@description:  write data to a csv file
*/

public class WriteToCsv {
	
	public static void writeToFile(String fileName, String str) {
		File results = new File(fileName);
		try {
			FileOutputStream fos = new FileOutputStream(results, true);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			osw.write(str + "\n");
			osw.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	* @Description: write array to excel file
	* @param data array
	* @param string  file name, must be xls, not xlsx   
	* @return a xls file   
	*/
	public void writeArrayExcel(double[][] data, String string) {
		int rowNum = data.length;
		int columnNum = data[0].length;
		try {
			FileWriter fw = new FileWriter(string);
			for (int i = 0; i < rowNum; i++) {
				for (int j = 0; j < columnNum; j++)
					fw.write(data[i][j]+ "\t"); // tab 间隔
				fw.write("\n"); // 换行
			}
			fw.close();
		}
		catch (IOException e){
			e.printStackTrace();
		}
	}
	
	/**
	* @Description: write an Arraylist to excel file
	* @param data arraylist
	* @param string  file name, must be xls, not xlsx   
	* @return a xls file   
	*/
	public void writeArrayExcel(ArrayList<double[]> data, String string) {
		int rowNum = data.size();
		int columnNum = data.get(0).length;
		try {
			FileWriter fw = new FileWriter(string);
			for (int i = 0; i < rowNum; i++) {
				for (int j = 0; j < columnNum; j++)
					fw.write(data.get(i)[j]+ "\t"); // tab 间隔
				fw.write("\n"); // 换行
			}
			fw.close();
		}
		catch (IOException e){
			e.printStackTrace();
		}
	}

}
