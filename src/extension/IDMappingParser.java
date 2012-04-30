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
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class IDMappingParser {
	
	ArrayList<String> uniprot_IDs = new ArrayList<String>();
	ArrayList<String> gene_names = new ArrayList<String>();
	ArrayList<String> gene_desc = new ArrayList<String>();
	
	public void parseIDMapping(String infile) {
		
		String line;
		StringTokenizer strtok;
		
		try {
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 br.readLine();
			 
			 while ((line = br.readLine()) != null) {
				 
				 strtok = new StringTokenizer(line, "\t");
				 uniprot_IDs.add(strtok.nextToken().trim());
				 strtok.nextToken();
				 gene_desc.add(strtok.nextToken().trim());
				 gene_names.add(strtok.nextToken().trim().toUpperCase());	 
			 }
			 br.close(); 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing species mapping.");
		}	
	}
	
public void writeMapping2File(String outfile) {
		
		try {
			 BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			 
			 for (int i=0; i<gene_names.size(); i++) {
				 
				 bw.write(gene_names.get(i) + "\n");
			 }
			 bw.flush();
			 bw.close(); 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing species mapping.");
		}	
	}
	
	public static void main(String[] args) {
		
		String mapdir = "/afs/cs/ra/share/eichner/pwm_procaryotes/MappingFiles/";
		String mapfile = "eColi_K12_MappingFile.txt";
		
		IDMappingParser mapping_parser = new IDMappingParser();
		mapping_parser.parseIDMapping(mapdir + mapfile);
		mapping_parser.writeMapping2File(mapdir + "parsed_genenames.txt");
	}
}


