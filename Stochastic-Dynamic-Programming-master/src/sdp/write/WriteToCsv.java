package sdp.write;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

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

}
