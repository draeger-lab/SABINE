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
import java.util.ArrayList;
import java.util.StringTokenizer;




/*
 * 
 * scales all files in given directory with svm_scale
 * 
 */
public class LibSVMFeatureScaler {
	
	
	/*
	 * removes "TF i vs. TF j : " - labeling of LibSVM-formatted files
	 */
	
	public void adjustSingleFeatureFile(String infile, String outfile) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(infile)));
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			String line = null;
			
			StringTokenizer strtok = null;
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line, ":");
				
				strtok.nextToken();
				
				strtok = new StringTokenizer(strtok.nextToken("x"));
				
				strtok.nextToken();
				
				bw.write(strtok.nextToken("x").trim() + "\n");
				
			}
			
			bw.flush();
			bw.close();
			br.close();

		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
		}
		
	}
	

	
	
	public void scaleFeatureFile(String orientationfile, String infile, String outfile) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		String line = null;
		
		StringTokenizer strtok = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(orientationfile)));
			
			ArrayList<ArrayList<String>> orientation_features = new ArrayList<ArrayList<String>>();
			
			boolean first = true;
			
			int num_features = 0;
			
			String[] split = null;
			
			
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				ArrayList<String> a_line = new ArrayList<String>();
				
				
			// count and parse features of first data point
				
				if(first) {
					
					a_line.add(strtok.nextToken());		// label
					
					while(strtok.hasMoreTokens()) {
						
						num_features++;
						
						split = strtok.nextToken().split(":");
						
						if(!split[0].equals(num_features + "")) {
							
							System.out.println("Error while parsing features of orientationfile. " + (num_features) + " expected. Aborting.");
							System.out.println("Line: " + line);
							System.exit(0);
							
						}
						
						a_line.add(split[1]);			// feature
						
						
					}
					
					first = false;
					
				}
				
				
			// parse features of another data point
				
				else {
					
					a_line.add(strtok.nextToken());		// label
					
					for(int i=0; i<num_features; i++) {
						
						split = strtok.nextToken().split(":");
						
						if(!split[0].equals((i+1) + "")) {
							
							System.out.println("Error while parsing features of orientationfile. " + (i+1) + " expected. Aborting.");
							System.out.println("Line: " + line);
							System.exit(0);
							
						}
						a_line.add(split[1]);			// feature
						
					}
					
				}
				
				orientation_features.add(a_line);
				
			}
			
			br.close();
			
			
			/*
			 * 
			 * fill double matrix with feature values + calculate min & max values
			 * 
			 */
			
			
			double[] max = new double[num_features];
			double[] min = new double[num_features];
			
			
			double[][] feature_matrix = new double[orientation_features.size()][num_features];
			
			for(int j=0; j<num_features; j++) {
				
				max[j] = Double.NEGATIVE_INFINITY;
				min[j] = Double.POSITIVE_INFINITY;
				
				for(int i=0; i<feature_matrix.length; i++) {
				
					feature_matrix[i][j] = Double.parseDouble(orientation_features.get(i).get(j+1));
					
					if(feature_matrix[i][j] < min[j]) min[j] = feature_matrix[i][j];
					if(feature_matrix[i][j] > max[j]) max[j] = feature_matrix[i][j];
					
					
				}
			}
			
			
			/*
			 * 
			 * calculate scaling parameters
			 * 
			 */
			
			double[] multiply 	= new double[num_features];
			double[] add		= new double[num_features];
			
			for(int j=0; j<num_features; j++) {
				
				add[j] = 0.0 - min[j];
				
				multiply[j] = 2.0 / ( max[j] - min[j] ); 
				
				if(max[j] == min[j]) {
					
					multiply[j] = 1.0 / max[j];
					
					if(max[j] == 0.0) {
						
						multiply[j] = 1.0;
						
					}
					
				}
			}
			
			
			/*
			 * 
			 * scale infile with these params
			 * 
			 */
			
			
			br = new BufferedReader(new FileReader(new File(infile)));
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				bw.write(strtok.nextToken());		// label
				
				
			// parse, scale and print features of a data point
					
				for(int i=0; i<num_features; i++) {
						
					split = strtok.nextToken().split(":");
					
					if(!split[0].equals((i+1) + "")) {
							
						System.out.println("Error while parsing features of infile. " + (i+1) + " expected. Aborting.");
						System.out.println("Line: " + line);
						System.exit(0);
							
					}
						
					bw.write( " " + (i+1) + ":" + ( ( Double.parseDouble(split[1]) + add[i] ) * multiply[i]  - 1.0));
						
				}
				
				bw.write("\n");
					
				
			}
			
			br.close();
			
			bw.flush();
			bw.close();
			
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while scaling feature file.");
		}
		
	}
}

