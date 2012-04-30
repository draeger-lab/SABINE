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

package optimization;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;


/*
 * 
 * constructs the reverse complement of a given PWM
 * 
 */

public class MatrixConverter {
	
	public String calculateConsensusLetter(float a, float c, float g, float t) {
		
		if(a > 0.6) return "A";
		if(c > 0.6) return "C";
		if(g > 0.6) return "G";
		if(t > 0.6) return "T";
		
		if(a + c > 0.8) return "M";
		if(a + g > 0.8) return "R";
		if(a + t > 0.8) return "W";
		if(c + g > 0.8) return "S";
		if(c + t > 0.8) return "Y";
		if(g + t > 0.8) return "K";
		
		return "N";
		
	}
	
	/*
	 * 
	 * gets a matrix in STAMP format and returns its reverse complement
	 * 
	 */
	public ArrayList<String> generateReverseComplement(ArrayList<String> matrix) {
		
	/*	
	  
	 	System.out.println("\nGenerating reverse complement.\n");
		
		System.out.println("\nMatrix:\n");
		
		for(int i=0; i<matrix.size(); i++) {
			System.out.println("Line : " + matrix.get(i));
		}
		
		System.out.println("\n");
	 
	 */
		
		
		StringTokenizer strtok = null;
		
		ArrayList<String> res = new ArrayList<String>();
		
		String a,c,g,t;
		
		String line = null;
		
		for(int i=matrix.size()-1; i>=0; i--) {
		
	//		System.out.println("Reversed : " + matrix.get(i));
			
			 strtok = new StringTokenizer(matrix.get(i));
			 
			 strtok.nextToken();
			 
			 a = strtok.nextToken();
			 c = strtok.nextToken();
			 g = strtok.nextToken();
			 t = strtok.nextToken();
			 
			 line = t + "\t" + g + "\t" + c + "\t" + a;
			 
			 res.add(line);
			 
		}
		
	//	System.out.println("Finished.");
		
		return res;
		
	}
	
	public void generateReverseComplements(String infile, String outfile, String format) {
		
		if(!format.equals("Matlign") && !format.equals("STAMP")) {
			
			System.out.println("Unknown Format.");
			System.out.println("Known Formats: STAMP , Matlign.");
			System.exit(0);
			
		}
		
		ArrayList<int[]> matrix = null;
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File("interleaved" + File.separator + infile)));
			bw = new BufferedWriter(new FileWriter(new File("interleaved" + File.separator + outfile)));
			
			String line = null;
			StringTokenizer strtok = null;
			
			
			
		// go over all matrices	
			
			while((line = br.readLine()) != null) {
				
				System.out.println(line);
				
				if(format.equals("Matlign")) {
					if(!line.startsWith(">")) {
						System.out.println("Parse Error. > expected.");
						System.exit(0);
					}
					else {
						bw.write(line + "_reversed\n");
					}
				}
				
				if(format.equals("STAMP")) {
					if(!line.startsWith("DE")) {
						System.out.println("Parse Error. DE expected.");
						System.exit(0);
					}
					else {
						bw.write(line + "_reversed\n");
					}
					
				}
					
			
			// parse current matrix
				
				matrix = new ArrayList<int[]>();
				
				
				if(format.equals("Matlign")) {
					
					while((line = br.readLine()) != null && line.length() > 0) {
						
						strtok = new StringTokenizer(line);
						
						matrix.add (
								
								new int[] {
										
										Integer.parseInt(strtok.nextToken()),
										Integer.parseInt(strtok.nextToken()),
										Integer.parseInt(strtok.nextToken()),
										Integer.parseInt(strtok.nextToken())
										
								}
								
						);
						
						
					}
				}
				
				if(format.equals("STAMP")) {
					
					while((line = br.readLine()) != null && !line.startsWith("XX")) {
						
						strtok = new StringTokenizer(line);
						
						strtok.nextToken();
						
						matrix.add (
								
								new int[] {
										
										Integer.parseInt(strtok.nextToken()),
										Integer.parseInt(strtok.nextToken()),
										Integer.parseInt(strtok.nextToken()),
										Integer.parseInt(strtok.nextToken())
										
								}
								
						);
						
					}
				}
				
				
			// reverse matrix and write it to output file
				
				if(format.equals("Matlign")) {
				
					for(int i=matrix.size()-1; i>=0; i--) {
						
						bw.write(matrix.get(i)[3] + "\t");
						bw.write(matrix.get(i)[2] + "\t");
						bw.write(matrix.get(i)[1] + "\t");
						bw.write(matrix.get(i)[0] + "\n");
					}
					
					bw.write("\n");
				}
				
				if(format.equals("STAMP")) {
					
					int line_count = 0;
					float a,c,g,t;
					float divisor = 100;
					
					for(int i=matrix.size()-1; i>=0; i--) {
						
						bw.write((line_count++)   + "\t");
						
						bw.write(matrix.get(i)[3] + "\t");
						bw.write(matrix.get(i)[2] + "\t");
						bw.write(matrix.get(i)[1] + "\t");
						bw.write(matrix.get(i)[0] + "\t");
						
						a = (float) matrix.get(i)[3] / divisor;
						c = (float) matrix.get(i)[2] / divisor;
						g = (float) matrix.get(i)[1] / divisor;
						t = (float) matrix.get(i)[0] / divisor;
						
						bw.write(calculateConsensusLetter(a, c, g, t) + "\n");
					
					}
					
					bw.write("XX\n");
					
				}
				
			}
			
			bw.flush();
			bw.close();
			br.close();
			
		}
		catch(FileNotFoundException fnfe) {
			System.out.println("File not found.");
		}
		catch(IOException ioe) {
			System.out.println("IO Exception occured.");
		}	
		
	}
	
	public static void main(String[] args) {
		
		MatrixConverter converter = new MatrixConverter();
		converter.generateReverseComplements("TestMatrix", "TestMatrix_reversed.out", "STAMP");
	
	}
	
}

