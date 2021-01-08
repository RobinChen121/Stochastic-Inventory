package sdp.write;

import java.io.File;
import java.io.FileInputStream;
import java.util.Arrays;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;



public class ReadExcel  
{  
	
	// hssf is for reading xls files
	public double[][] readExcelXLSX(String fileAddress, int rowStartIndex) {
		try {
			File file = new File(fileAddress); // creating a new file instance
			FileInputStream fis = new FileInputStream(file); // obtaining bytes from the file
			
			// creating Workbook instance that refers to .xlsx file
			XSSFWorkbook wb = new XSSFWorkbook(fis);
			XSSFSheet st = wb.getSheetAt(0); // creating a Sheet object to retrieve object
			
			double[][] paraTable = new double[st.getLastRowNum()- rowStartIndex + 1][];
			// starting row number
			for (int rowIndex = rowStartIndex; rowIndex <= st.getLastRowNum(); rowIndex++) {
				 XSSFRow row = st.getRow(rowIndex);
				 paraTable[rowIndex- rowStartIndex] = new double[row.getLastCellNum()];
				 for (int columnIndex = 0; columnIndex <= row.getLastCellNum() - 1; columnIndex++) {
					 Cell cell = row.getCell(columnIndex);
					 paraTable[rowIndex - rowStartIndex][columnIndex] = cell.getNumericCellValue();
				 }
			}
			
			wb.close();
			return paraTable;
			
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	
	public static void main(String[] args) {
		double[][] test = new ReadExcel().readExcelXLSX("Numerical experiments-settings.xlsx", 2);
		System.out.println(Arrays.deepToString(test));
	}
}
