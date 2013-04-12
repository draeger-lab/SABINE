package extension;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class BasicTools {

	/**
	 * 
	 * @param fileName
	 * @param upperCase
	 * @return
	 */
	public static List<String> readFile2List(String fileName, boolean upperCase) {
		
		List<String> fileContent = null;
		
		try {
			fileContent = readStream2List(new FileInputStream(new File(fileName)), upperCase);
		
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return fileContent;
	}
	
	/**
	 * 
	 * @param stream
	 * @param upperCase
	 * @return
	 */
	public static List<String> readStream2List(InputStream stream, boolean upperCase) {
		
		List<String> fileContent = new ArrayList<String>();
		
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(stream));
			
			String line;
			while ((line = br.readLine()) != null) {
				if (upperCase) {
					fileContent.add(line.trim().toUpperCase());
				} else {
					fileContent.add(line.trim());
				}
			}
			br.close();
			
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		
		return(fileContent);
	}
	
	

	public static void writeList2File(List<String> list, String outfile) {
		String[] array = list.toArray(new String[] {});
		writeArray2File(array, outfile);
	}

	public static void writeArray2File(String[] array, String outfile) {
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		
			for (String line: array) {
				bw.write(line + "\n");
			}
			bw.flush();
			bw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
