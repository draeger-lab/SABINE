package cv;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

public class CVHelper {

	/**
	 * 
	 * For a given gram matrix check if this is a classification or regression
	 * problem
	 * 
	 * @param f
	 * @return true, if regression problem
	 * @throws IOException
	 */
	@SuppressWarnings("deprecation")
	public static boolean isARegressionProblem(File f) throws IOException {
		ArrayList<String> al = new ArrayList<String>();
		try {
			DataInputStream read = new DataInputStream(new BufferedInputStream(
					new FileInputStream(f)));
			String r;
			while ((r = read.readLine()) != null) {
				int pos = r.indexOf(" ");
				if(pos < 0){
					  pos = r.indexOf("\t");
				}
				String sub = r.substring(0, pos);
				// System.out.println(sub);
				if (!al.contains(sub)) {
					al.add(sub);
					if (al.size() > CVGlobalSettings.regression_threshold) {
						return true;
					}
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return false;
	}

	/**
	 * 
	 * simple check if a libsvm file is corrupt
	 * 
	 * @param f
	 * @return
	 * @throws IOException
	 */
	@SuppressWarnings("deprecation")
	public static boolean isEmptyFile(File f) throws IOException {
		try {
			DataInputStream read = new DataInputStream(new BufferedInputStream(
					new FileInputStream(f)));
			int max_line_length = 0;
			String line = "";
			while ((line = read.readLine()) != null) {
				if(line.length() > max_line_length){
					max_line_length = line.length();
				}
			}
			/**
			 * return false, if line length never exceeds 1
			 * this should detect most corrupt input files
			 */
			if(max_line_length < 2){
				return true;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return false;
	}
}
