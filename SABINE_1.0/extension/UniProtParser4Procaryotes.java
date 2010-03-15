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
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class UniProtParser4Procaryotes {

	ArrayList<ArrayList<String>> uniprot_IDs = new ArrayList<ArrayList<String>>();
	ArrayList<String> sequences = new ArrayList<String>();
	ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
	
	ArrayList<String> all_uniprot_IDs = new ArrayList<String>();
	ArrayList<Integer[]> id_map = new ArrayList<Integer[]>();

	public void parseSequencesAndDomains(String infile) {
		
		String line, curr_seq;
		StringTokenizer strtok;
		ArrayList<String> curr_IDs = new ArrayList<String>();
		ArrayList<String> curr_domains = new ArrayList<String>();
		boolean id_parsed, seq_parsed;
		id_parsed = seq_parsed = false;
		
		int line_cnt = 0;
		
		try { 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 
			 while ((line = br.readLine()) != null) {
				 
				 line_cnt++;
				 
				 // parse UniProt ID
				 while (line.startsWith("AC")) {
					 strtok = new StringTokenizer(line.substring(5), ";");
					 
					 while (strtok.hasMoreTokens()) {
						 curr_IDs.add(strtok.nextToken());
					 }
					 line = br.readLine();
					 line_cnt++;
					 id_parsed = true;
				 }
				 
				 // parse DNA binding domain
				 while (line.startsWith("FT   DNA_BIND")) {
					 strtok = new StringTokenizer(line.substring(14));
					 curr_domains.add(strtok.nextToken() + "  " + strtok.nextToken());
					 line = br.readLine();
					 line_cnt++;
				 }
				 
				 // parse sequence
				 if (line.startsWith("SQ")) {
					 curr_seq = "";
					 seq_parsed = true;
					 while (! (line=br.readLine()).startsWith("//")) {
						 curr_seq += line.trim().replace(" ", "");
						 line_cnt++;
					 }
					 line_cnt++;
					 sequences.add(curr_seq);
				 }
				 
				 if (line.startsWith("//")) {
					 if (id_parsed && seq_parsed) {
						 id_parsed = seq_parsed = false;
						 
						 uniprot_IDs.add(curr_IDs);
						 curr_IDs = new ArrayList<String>();
						 
						 domains.add(curr_domains);
						 curr_domains = new ArrayList<String>();
					 }
					 else {
						 System.out.println("Error occurred while parsing sequences and domains from UniProt.");
						 System.exit(0);
					 }
				 }
			 }
			for (int i=0; i<uniprot_IDs.size(); i++) {
				for (int j=0; j<uniprot_IDs.get(i).size(); j++) {
						
					all_uniprot_IDs.add(uniprot_IDs.get(i).get(j));
					id_map.add(new Integer[] {i,j});
				}
			}
			//System.out.println("  " + uniprot_IDs.size() + " UniProt IDs, domains and sequences parsed.");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing sequences and domains from UniProt.");
		}
	}
	
	
	public static void main(String[] args) {
		
		//String path = "/afs/cs/ra/share/eichner/pwm_procaryotes/Uniprot_Annotations/";
		String path = "S:/share/eichner/pwm_procaryotes/Uniprot_Annotations/";
		String infile = path + "EscherichiaColi(strain K12)_ECOLI.dat";
		
		UniProtParser4Procaryotes uniprot_parser = new UniProtParser4Procaryotes();
		uniprot_parser.parseSequencesAndDomains(infile);
	}
}

