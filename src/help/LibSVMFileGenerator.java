/*
    SABINE predicts binding specificities of transcription factors.
    Copyright (C) 2009 ZBIT, University of Tübingen, Johannes Eichner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package help;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.StringTokenizer;

public class LibSVMFileGenerator {
	
	
	public void generateUnlabeledLibSVMFile(String featurefile, String outfile) {
		
		
		DecimalFormat fmt = new DecimalFormat();
		
		fmt.setMaximumFractionDigits(16);
		fmt.setMinimumFractionDigits(16);
			
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
		symbs.setDecimalSeparator('.');
		fmt.setDecimalFormatSymbols(symbs);
		
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		try {
			
			
			br = new BufferedReader(new FileReader(new File(featurefile)));
			bw  = new BufferedWriter(new FileWriter(new File(outfile)));
			
			String line = null;
			
			StringTokenizer strtok = null;
			
			
		// go over all labels
			 
			while((line = br.readLine()) != null) {
				
				
				strtok = new StringTokenizer(line);
				
				strtok.nextToken(":");
				
				
				
			// write label and first feature	
				
				bw.write(fmt.format(0.0) + " 1:");
				bw.write(fmt.format(Double.parseDouble(strtok.nextToken(":").substring(1))).replace(",", "") + "\n");
				
				
			}
			
			
			br.close();
			bw.flush();
			bw.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while generating input file for libsvm.");			
		}
	
	}
	
	
public void generateLabeledLibSVMFile(String featurefile, String labelfile, String outfile) {
		
		
		DecimalFormat fmt = new DecimalFormat();
		
		fmt.setMaximumFractionDigits(16);
		fmt.setMinimumFractionDigits(16);
			
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
		symbs.setDecimalSeparator('.');
		fmt.setDecimalFormatSymbols(symbs);
		
		
		BufferedReader br1 = null;
		BufferedReader br2 = null;
		BufferedWriter bw = null;
		
		try {
			
			
			br1 = new BufferedReader(new FileReader(new File(featurefile)));
			br2 = new BufferedReader(new FileReader(new File(labelfile)));
			bw  = new BufferedWriter(new FileWriter(new File(outfile)));
			
			String line1 = null;
			String line2 = null;
			
			StringTokenizer strtok1 = null;
			StringTokenizer strtok2 = null;
			
			
		// go over all labels
			 
			while((line1 = br1.readLine()) != null) {
				
				if ((line2 = br2.readLine()) == null) {
					System.out.println("Parse Error. Feature file and label file have different sizes.");
					System.exit(0);
				}
				
				strtok1 = new StringTokenizer(line1);
				strtok2 = new StringTokenizer(line2);
				
				strtok1.nextToken(":");
				strtok2.nextToken(":");
				
				
			// write label and first feature	
				
				bw.write(fmt.format(Double.parseDouble(strtok2.nextToken(":").substring(1))).replace(",", "") + " 1:");
				bw.write(fmt.format(Double.parseDouble(strtok1.nextToken(":").substring(1))).replace(",", "") + "\n");
				
				
			}
			
			
			br1.close();
			br2.close();
			bw.flush();
			bw.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while generating input file for libsvm.");			
		}
	
	}
	
	
	
	public void addFeature(String svmfile, String featurefile) {
		
		BufferedReader br1 = null;
		BufferedReader br2 = null;
		BufferedWriter bw = null;
		
		DecimalFormat fmt = new DecimalFormat();
		
		fmt.setMaximumFractionDigits(16);
		fmt.setMinimumFractionDigits(16);
			
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
		symbs.setDecimalSeparator('.');
		fmt.setDecimalFormatSymbols(symbs);
		
		
		boolean first = true;
		
		String tmpdir = svmfile.substring(0,svmfile.lastIndexOf("/")+1); 
		
		try {
			
			br1 = new BufferedReader(new FileReader(new File(svmfile)));
			br2 = new BufferedReader(new FileReader(new File(featurefile)));
			bw  = new BufferedWriter(new FileWriter(new File(tmpdir + "tmp.out")));
			
			String line1 = null;
			String line2 = null; 
			
			StringTokenizer strtok2 = null;
			
			int featureIndex = 0;
			
		// go over all entries in libsvm-file
			 
			while((line1 = br1.readLine()) != null) {
				
				line2 = br2.readLine();
				
				if(line2 == null) {
					System.out.println("Error. More labels than features. Aborting.");
					System.exit(0);	
				}
				
				//System.out.println("Line: " + line1);
				
				
			// count number of old features	
				
				if(first) {
					
					first = false;
					StringTokenizer strtok1 = new StringTokenizer(line1);
					String token = null;
					
				// goto last feature	
					
					while(strtok1.hasMoreTokens()) {
						token = strtok1.nextToken();
					}
					
					
					String[] split = token.split(":");
					//System.out.println(split[0]);
					
					if(split.length != 2) {
						System.out.println("Error while counting number of features. Aborting.");
						System.exit(0);
					}
				
				// init feature index for new feature
					
					featureIndex = Integer.parseInt(split[0]) + 1;
					
				}
				
				strtok2 = new StringTokenizer(line2);
				strtok2.nextToken(":");
				
				
			// write label & old features and add new feature	
				
				bw.write(line1 + " " + featureIndex + ":");
				
				bw.write(fmt.format(Double.parseDouble(strtok2.nextToken(":").substring(1))).replace(",", "") + "\n");
				
				
			}
			
			if((line2 = br2.readLine()) != null) {
				System.out.println("Error. More features than labels. Aborting.");
				System.exit(0);	
			}
			
			
			br1.close();
			br2.close();
			bw.flush();
			bw.close();
			
			
			// overwrite old svmfile
			
			try {
				Process proc = Runtime.getRuntime().exec("mv " + tmpdir + "tmp.out " + svmfile);
				proc.waitFor();
				proc.destroy();
				
			} catch (InterruptedException e) {
				e.printStackTrace();
			}

			
			if (! (new File(svmfile)).exists())
				System.out.print("SVM-file was not found.");
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while adding features to input file for libsvm.");			
		}
		
	}
	
	
}

