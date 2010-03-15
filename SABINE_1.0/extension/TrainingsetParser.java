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

package extension;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

import model.LabelFileGenerator;

public class TrainingsetParser {

  public ArrayList<String> tf_names = new ArrayList<String>();
  ArrayList<String> sequences1 = new ArrayList<String>();
  ArrayList<String> sequences2 = new ArrayList<String>();
  ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
  
  ArrayList<String> species = new ArrayList<String>();
  public ArrayList<ArrayList<String>> matrices = new ArrayList<ArrayList<String>>();
  
  public void parseFactors(String infile) {
	  
	  String line, curr_seq;
	  StringTokenizer strtok;
	  ArrayList<String> curr_domains;
	  
	  try {
		  
		  BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
		  
		  while ((line = br.readLine()) != null) {
			  
			  /*
			   * parse factor name
			   */
			  
			  if (line.startsWith("NA  ")) {
				  tf_names.add(line.substring(4).trim());
			  }
			  
			  /*
			   *  parse DNA binding domains
			   */
			  
			  if (line.startsWith("DO  ")) {
				  
				  curr_domains = new ArrayList<String>();
				  while (line.startsWith("DO  ")) {
					
					  strtok = new StringTokenizer(line.substring(4));
					  
					  strtok.nextToken();
					  strtok.nextToken();
					  
					  curr_domains.add(strtok.nextToken().trim() + "\t" + strtok.nextToken().trim());
					  
					  line = br.readLine();
				  }
				  domains.add(curr_domains);
			  }
				  
			  /*
			   *  parse sequences
			   */
			  
			  if (line.startsWith("S1  ")) {
				  curr_seq = line.substring(4).trim();
				  while((line = br.readLine()).startsWith("S1  ")) {
					  curr_seq += line.substring(4).trim();
				  }
				  sequences1.add(curr_seq);
			  }
			  if (line.startsWith("S2  ")) {
				  curr_seq = line.substring(4).trim();
				  while((line = br.readLine()).startsWith("S2  ")) {
					  curr_seq += line.substring(4).trim();
				  }
				  sequences2.add(curr_seq);
			  }
			  
			  if (line.startsWith("//") && sequences1.size() < sequences2.size())
				  sequences1.add(null);
			  if (line.startsWith("//") && sequences1.size() > sequences2.size())
				  sequences2.add(null);
		  }
	  }
	  catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing training factors.");
	  }
  }
  
  public void parseSpecies(String infile) {
	
	  String line;
	  
	  try {
		BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
		  
		while ((line = br.readLine()) != null) {
			species.add(line.trim());
		}
	  }	
	  catch(IOException ioe) {
		System.out.println(ioe.getMessage());
		System.out.println("IOException occurred while parsing species of training factors.");
  	  }
  }
  
  public void parseMatrices(String infile) {
	  
	  LabelFileGenerator matrix_parser = new LabelFileGenerator();
	  matrices = matrix_parser.getAllFBPs(infile);  
  }
  
  public void parseAll(int class_id) {
	  
	  parseFactors("trainingsets/trainingset_class" + class_id + ".rawdata");
	  parseSpecies("trainingsets/trainingset_class" + class_id + ".species");
	  parseMatrices("trainingsets/trainingset_class" + class_id + ".fbps");
  }
  
  
	public void writeInputFile(BufferedWriter bw, int[] relevant_factors, int class_id) {
		
		String curr_seq;
		int SEQLINELENGTH = 60;
		int i;
		
		try {
			
			for (int x=0; x<relevant_factors.length; x++) {
				
				i = relevant_factors[x];
				
				/*
				 *    write name and species
				 */
				
				bw.write("NA  " + tf_names.get(i) + "\n");
				bw.write("XX\n" +   
						 "SP  " + species.get(i) + "\n" + 
						 "XX\n");
				
				/*
				 *    write TRANSFAC superclass
				 */
				
				bw.write("CL  " + class_id + ".1.1.1.1.\n" + 
						 "XX\n");
				
				/*
				 *    write sequence 
				 */
				
				curr_seq = sequences1.get(i);
				
				if (curr_seq != null) {
					for(int j=0; j<(curr_seq.length()/SEQLINELENGTH); j++) {
					
						bw.write("S1  "); 
						bw.write(curr_seq.toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
						bw.write("\n");
					}
					
					if(curr_seq.length()-(curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
						
						bw.write("S1  "); 
						bw.write(curr_seq.toUpperCase(), (curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH, curr_seq.length()-(curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");
					}
				
					bw.write("XX\n");
				}
				
				curr_seq = sequences2.get(i);
				
				if (curr_seq != null) {
					for(int j=0; j<(curr_seq.length()/SEQLINELENGTH); j++) {
					
						bw.write("S2  "); 
						bw.write(curr_seq.toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
						bw.write("\n");
					}
					
					if(curr_seq.length()-(curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
						
						bw.write("S2  "); 
						bw.write(curr_seq.toUpperCase(), (curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH, curr_seq.length()-(curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");
					}
				
					bw.write("XX\n");
				}
				
				/*
				 *    write domains 
				 */
				
				for (int j=0; j<domains.get(i).size(); j++) {
					bw.write("FT  " + "domainID" + "\t" + domains.get(i).get(j) + "\n");
				}
				
				bw.write("XX\n" + 
						 "//\n" +
						 "XX\n");	
			}
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing input file for SABINE.");
		}
	}
	
	
  
  public static void main(String[] args) {
	  
	  TrainingsetParser train_parser = new TrainingsetParser();
	  train_parser.parseAll(0);
  }

}

