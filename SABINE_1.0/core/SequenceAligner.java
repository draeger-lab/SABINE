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

package core;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.StringTokenizer;

import org.biojava.bio.BioException;
import org.biojava.bio.alignment.NeedlemanWunsch;
import org.biojava.bio.alignment.SequenceAlignment;
import org.biojava.bio.alignment.SmithWaterman;
import org.biojava.bio.alignment.SubstitutionMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;

 
/*
 * Substitution Matrix Files available at ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/
 */

public class SequenceAligner {
 
	
	/*
	 * 
	 * fields
	 * 
	 */
	
	FiniteAlphabet alphabet = null;
    
	SubstitutionMatrix matrix = null;
	
	SequenceAlignment alignment = null;
	
	StringTokenizer strtok = null;
	
	
	/*
	 * 
	 * constructors
	 * 
	 */
	
	public SequenceAligner(String matrixfile, String type) throws IOException, BioException {
		
		alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
		
		matrix = new SubstitutionMatrix(alphabet, new File("substitutionMatrices" + File.separator + matrixfile));
		
		if(type.equals("NW")) {
			alignment = new NeedlemanWunsch( 
			        0, 		// match
			        3,		// replace
			        11,  	// insert
			        11,		// delete
			        1,  	// gapExtend
			        matrix 	// SubstitutionMatrix
			      );
			      
		}
		else {
			alignment = new SmithWaterman(
			        -1,     // match
			        3,      // replace 
			        11,      // insert
			        11,      // delete
			        1,      // gapExtend
			        matrix  // SubstitutionMatrix
			      );
		}
	}
	
	
	public SequenceAligner(String matrixfile, int ins, int del, int gap_ext) throws IOException, BioException {
		
		alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
		
		matrix = new SubstitutionMatrix(alphabet, new File("substitutionMatrices" + File.separator + matrixfile));
		
		
		alignment = new NeedlemanWunsch( 
			        0, 			// match
			        3,			// replace
			        ins,  		// insert
			        del,		// delete
			        gap_ext,  	// gapExtend
			        matrix 		// SubstitutionMatrix
		);
			      
		
		
	}
	
	
	
	public SequenceAligner(int match, int mismatch, String type) throws IOException, BioException {
		
		alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
		
		matrix = new SubstitutionMatrix(alphabet, match, mismatch);
		
		if(type.equals("NW")) {
			alignment = new NeedlemanWunsch( 
			        -1, 	// match
			        1,		// replace
			        11,  	// insert
			        11,		// delete
			        1,  	// gapExtend
			        matrix 					// SubstitutionMatrix
			      );
			      
		}
		else {
			alignment = new SmithWaterman(
			        0,      // match
			        1,      // replace 
			        11,      // insert
			        11,      // delete
			        1,      // gapExtend
			        matrix  // SubstitutionMatrix
			      );
		}
	
	}
	
	
	public SequenceAligner(String dna, String matrixfile, String type) throws IOException, BioException {
		
		if(!dna.equals("DNA")) {
			
			System.out.println("First Argument must be \"DNA\". Aborting.");
			
			System.exit(0);
			
		}
		
		
		alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
		
		matrix = new SubstitutionMatrix(alphabet, new File("substitutionMatrices" + File.separator + matrixfile));
		
		
		if(type.equals("NW")) {
			alignment = new NeedlemanWunsch( 
			        0, 		// match
			        0,		// replace
			        10,  	// insert
			        10,		// delete
			        1,  	// gapExtend
			        matrix 	// SubstitutionMatrix
			      );
			      
		}
		
		if(type.equals("NW_UNGAPPED")) {
			alignment = new NeedlemanWunsch( 
			        0, 						// match
			        0,						// replace
			        Integer.MAX_VALUE,  	// insert
			        Integer.MAX_VALUE,		// delete
			        Integer.MAX_VALUE,  	// gapExtend
			        matrix 					// SubstitutionMatrix
			      );
			      
		}
		
		if(type.equals("SW")) {
			alignment = new SmithWaterman(
			        0,     // match
			        0,      // replace 
			        10,      // insert
			        10,      // delete
			        1,      // gapExtend
			        matrix  // SubstitutionMatrix
			      );
		}
	}
	
	
	/*
	 * 
	 * methods
	 * 
	 */
	
	public double getSequenceSimilarity(String seq1, String seq2) {
	
		double similarity = 0.0;
		
		String upperseq = null;
		String dashes = null;
		
		StringTokenizer strtok1 = null;
		StringTokenizer strtok2 = null;
		
		try {
			
			Sequence query  = ProteinTools.createProteinSequence(seq1, "query");
		    Sequence target = ProteinTools.createProteinSequence(seq2, "target");
		      
		    alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString(),"\n");
		
		}
		catch (Exception e) {
			System.out.println("Exception occurred while trying to align sequences.");
		}
		

		DecimalFormat fmt = new DecimalFormat();
		
		fmt.setMaximumFractionDigits(2);
		fmt.setMinimumFractionDigits(2);
			
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
		symbs.setDecimalSeparator('.');
		fmt.setDecimalFormatSymbols(symbs);
		
		String line = null;
		
		
		
	// goto start of alignment	
		
		for(int i=0; i<6; i++) line = strtok1.nextToken();
		
		
		if(!line.startsWith("Query:")) {
			System.out.println("Error while parsing alignment string. Aborting.");
			System.exit(0);
		}
		
		int count_all 	= 0;
		int count_match = 0;
		
		
		
	// go over all lines of alignment	
		
		while(line.startsWith("Query:")) {
			
			strtok2 = new StringTokenizer(line.replace("Query:", ""));
			strtok2.nextToken();
			
			upperseq = strtok2.nextToken();
			
							
			line = strtok1.nextToken(); // goto next line
			
			
			dashes = line.substring(line.length()- upperseq.length() - 1,line.length()-1); // trim whitespaces before |'s
			
			
			dashes = dashes.substring(0,upperseq.length()); // trim whitespace after |'s
			
			
			
			for(int i=0; i<dashes.length(); i++) {
				
				count_all++;
				
				if(dashes.charAt(i) == '|') 
					count_match++;
			}
			
		if(strtok1.hasMoreTokens())	line = strtok1.nextToken(); // goto next line
		if(strtok1.hasMoreTokens())	line = strtok1.nextToken(); // goto next line
		if(strtok1.hasMoreTokens())	line = strtok1.nextToken(); // goto next line
		
			
		}
		
		similarity = (double) count_match / (double) count_all;
		
		
		return similarity;
	}
	
	
	/*
	 * 
	 * calculates a sequence identity which ignores mismatches whose SM-Scores are above a given threshold
	 * 
	 */
	
	public double getSMBasedIdentity(String seq1, String seq2, double threshold) {
	
		
		double similarity = 0.0;
		

		try {
			
			/*
			 * 
			 * align sequences and calculate score
			 * 
			 */
						
			Sequence query  = ProteinTools.createProteinSequence(seq1, "query");
		    Sequence target = ProteinTools.createProteinSequence(seq2, "target");
		      
		    alignment.pairwiseAlignment(query, target);
		    
		    
		    StringTokenizer strtok  = new StringTokenizer(alignment.getAlignmentString(), "\n");
		    StringTokenizer strtok2 = null;
		    
		    
		    
		    String line = null;
		    
		    String upperseq = null;
		    String lowerseq = null;
		    String dashes   = null;
			
		    		    
		    
			// goto start of alignment	
			for(int i=0; i<6; i++) line = strtok.nextToken();
		    
			if(!line.startsWith("Query:")) {
				System.out.println("Error while parsing alignment string. Aborting.");
				System.exit(0);
			}
				
			int count_all 	= 0;
			int count_match = 0;
		
			
			while(line.startsWith("Query:")) {
				
				strtok2 = new StringTokenizer(line.replace("Query:", ""));
				strtok2.nextToken();
				
				upperseq = strtok2.nextToken();
				
								
				line = strtok.nextToken(); // goto next line
				
				
				dashes = line.substring(line.length()- upperseq.length() - 1,line.length()-1); // trim whitespaces before |'s
				
				
				dashes = dashes.substring(0,upperseq.length()); // trim whitespace after |'s
				
				
				line = strtok.nextToken(); // goto next line
				
				strtok2 = new StringTokenizer(line.replace("Target:", ""));
				strtok2.nextToken();
				
				lowerseq = strtok2.nextToken();
				
			//	System.out.println(upperseq + "\n" + dashes + "\n" + lowerseq + "\n");
				
				if(lowerseq.length() != dashes.length() || upperseq.length() != dashes.length()) {
					
					System.out.println("Error. Different lengthes of upper/lower sequence and dashes.\n");
					System.out.println("Upperseq : " + upperseq);
					System.out.println("Dashseq  : " + dashes);
					System.out.println("Lowerseq : " + lowerseq);
					System.exit(0);
						
				}
				
				Symbol upperSymbol = null;
				Symbol lowerSymbol = null;
				
				for(int i=0; i<dashes.length(); i++) {
					
					count_all++;
					
					boolean gap = false;
					
					if(upperseq.charAt(i) == '~' || lowerseq.charAt(i) == '~') gap = true;
					
					if(upperseq.charAt(i) == '-' || lowerseq.charAt(i) == '-') gap = true;
					
					if(!gap) {
						upperSymbol = getProteinSymbol(upperseq.charAt(i));
						lowerSymbol = getProteinSymbol(lowerseq.charAt(i));
					}
					
					
					if(	!gap && ( 

						dashes.charAt(i) == '|' 
						
							||
							
						matrix.getValueAt(upperSymbol, lowerSymbol) >= threshold ) ) 
						
						
						count_match++;
				}
				
			
				if(strtok.hasMoreTokens())	line = strtok.nextToken(); // goto next line
				if(strtok.hasMoreTokens())	line = strtok.nextToken(); // goto next line
			
			} 
			
			similarity = (double) count_match / (double) count_all;
		
		}
		catch(BioException be) {
			System.out.println(be.getMessage());
			System.out.println("BioException occurred while comparing sequences.");
		} 
		catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.out.println("Exception occurred while comparing sequences.");
			System.exit(0);
		}
		
		return similarity;
	}
	
	
	/*
	 * 
	 * calculates the substitution matrix based sequence similarity of two sequences
	 * 
	 */
	
	public double getSMBasedSimilarity(String seq1, String seq2) {
		
		
		double similarity = 0.0;
		
		double self_score1 = 0.0;
		double self_score2 = 0.0;
		
		StringTokenizer strtok1 = null;
		

		try {
			
			/*
			 * 
			 * align sequences and calculate score
			 * 
			 */
			
			Sequence query  = ProteinTools.createProteinSequence(seq1, "query");
		    Sequence target = ProteinTools.createProteinSequence(seq2, "target");
		      
		    alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());

			while(!strtok1.nextToken().equals("Score:"));
			
			similarity = Double.parseDouble(strtok1.nextToken());
			
			
			/*
			 * 
			 * normalize similarity score
			 * 
			 */
			
			
		// compare seq1 vs. itself	
			
			query  = ProteinTools.createProteinSequence(seq1, "query");
		    target = ProteinTools.createProteinSequence(seq1, "target");
		      
		    alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			self_score1 = Double.parseDouble(strtok1.nextToken());
			
			
		// compare seq2 vs. itself	
			
			query  = ProteinTools.createProteinSequence(seq2, "query");
		    target = ProteinTools.createProteinSequence(seq2, "target");
		      
		   alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			self_score2 = Double.parseDouble(strtok1.nextToken());
			
			
		// divide sequence similarity by atrithmetic mean of self scores
			
			similarity = 2 * similarity / (self_score1 + self_score2);
			
			
		}
		catch(BioException be) {
			System.out.println(be.getMessage());
		} 
		catch (Exception e) {
			e.getMessage();
			e.printStackTrace();
		}
		
		
		return similarity;
	}
	
	
	/*
	 *  returns the alignment score and the matching region of the query
	 */
	
	public String[] getMatchingRegionAndScore(String seq1, String seq2) {
		
		double similarity = 0.0;
		
		double self_score1 = 0.0;
		double self_score2 = 0.0;
		
		StringTokenizer strtok1 = null;
		
		String[] scoreAndRegion = new String[3];

		try {
			
			/*
			 * 
			 * align sequences and calculate score
			 * 
			 */
			
			Sequence query  = ProteinTools.createProteinSequence(seq1, "query");
		    Sequence target = ProteinTools.createProteinSequence(seq2, "target");
		      
		    alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			similarity = Double.parseDouble(strtok1.nextToken());
			
			
			/*
			 *  parse matching region of query
			 */
			
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			// initialize start and end position of matching region
			while(!strtok1.nextToken().equals("Query:"));   // skip output of query length
			while(!strtok1.nextToken().equals("Query:"));	// go to alignment (query sequence)
			String start_pos = strtok1.nextToken().trim(); 	// start position
			strtok1.nextToken();							// alignment
			String end_pos = strtok1.nextToken().trim();	// end position
			
			// find end position of matching region		
			while(strtok1.hasMoreTokens()) {
				if (strtok1.nextToken().equals("Query:")) {
					strtok1.nextToken();					// start position
					strtok1.nextToken();					// alignment
					end_pos = strtok1.nextToken().trim(); 	// end position
				}
			}
			
			/*
			 * 
			 * normalize similarity score
			 * 
			 */
			
			
		// compare seq1 vs. itself	
			
			query  = ProteinTools.createProteinSequence(seq1, "query");
		    target = ProteinTools.createProteinSequence(seq1, "target");
		      
		    alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			self_score1 = Double.parseDouble(strtok1.nextToken());
			
			
		// compare seq2 vs. itself	
			
			query  = ProteinTools.createProteinSequence(seq2, "query");
		    target = ProteinTools.createProteinSequence(seq2, "target");
		      
		   alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			self_score2 = Double.parseDouble(strtok1.nextToken());
			
			
		// divide sequence similarity by atrithmetic mean of self scores
			
			similarity = 2 * similarity / (self_score1 + self_score2);
			
			scoreAndRegion = new String[] { similarity + "" , start_pos , end_pos };
			/*
			System.out.println("Alignment score: " + similarity);
			System.out.println("Start position:  " + start_pos);
			System.out.println("End position:    " + end_pos);
			*/
			
		}
		catch(BioException be) {
			System.out.println(be.getMessage());
		} 
		catch (Exception e) {
			e.getMessage();
			e.printStackTrace();
		}
		
		
		
		
		return scoreAndRegion;
	}
	
	
	/*
	 * 
	 * calculates the substitution matrix based sequence similarity of two sequences
	 * 
	 */
	
	public double getSMBasedSecondaryStructureSimilarity(String seq1, String seq2) {
		
		
		double similarity = 0.0;
		
		double self_score1 = 0.0;
		double self_score2 = 0.0;
		
		StringTokenizer strtok1 = null;
		

		try {
			
			/*
			 * 
			 * align sequences and calculate score
			 * 
			 */
			
			Sequence query  = DNATools.createDNASequence(seq1, "query");
		    Sequence target = DNATools.createDNASequence(seq2, "target");
		      
		    alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			similarity = Double.parseDouble(strtok1.nextToken());
			
			
			/*
			 * 
			 * normalize similarity score
			 * 
			 */
			
			
		// compare seq1 vs. itself	
			
			query  = DNATools.createDNASequence(seq1, "query");
		    target = DNATools.createDNASequence(seq1, "target");
		      
		    alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			self_score1 = Double.parseDouble(strtok1.nextToken());
			
			
		// compare seq2 vs. itself	
			
			query  = DNATools.createDNASequence(seq2, "query");
		    target = DNATools.createDNASequence(seq2, "target");
		      
		   alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			self_score2 = Double.parseDouble(strtok1.nextToken());
			
			
		// divide sequence similarity by geometric mean of self scores
			
			similarity = similarity / Math.sqrt(self_score1 * self_score2);
			
			
		}
		catch(BioException be) {
			System.out.println(be.getMessage());
		} 
		catch (Exception e) {
			e.getMessage();
			e.printStackTrace();
		}
		
		
		return similarity;
	}
	
	
	
	public static Symbol getProteinSymbol(char c) {
		
		Symbol res = null;
		
		if (c == 'X') { c = 'G'; }
		
		switch(c) {
		
			case 'A' : { res = ProteinTools.a(); break; }
			case 'R' : { res = ProteinTools.r(); break; }
			case 'N' : { res = ProteinTools.n(); break; }
			case 'D' : { res = ProteinTools.d(); break; }
			case 'C' : { res = ProteinTools.c(); break; }
			case 'Q' : { res = ProteinTools.q(); break; }
			case 'E' : { res = ProteinTools.e(); break; }
			case 'G' : { res = ProteinTools.g(); break; }
			case 'H' : { res = ProteinTools.h(); break; }
			case 'I' : { res = ProteinTools.i(); break; }
			case 'L' : { res = ProteinTools.l(); break; }
			case 'K' : { res = ProteinTools.k(); break; }
			case 'M' : { res = ProteinTools.m(); break; }
			case 'F' : { res = ProteinTools.f(); break; }
			case 'P' : { res = ProteinTools.p(); break; }
			case 'S' : { res = ProteinTools.s(); break; }
			case 'T' : { res = ProteinTools.t(); break; }
			case 'W' : { res = ProteinTools.w(); break; }
			case 'Y' : { res = ProteinTools.y(); break; }
			case 'V' : { res = ProteinTools.v(); break; }
			
			default  : { System.out.println("Illegal Symbol: " + c + ". Aborting. "); System.exit(0); break; }
			
		}
		
		return res;
		
	}
	

	public double getLAKernelScore(String seq1, String seq2, String matrix) {
		
		
		if(matrix == null) matrix = "blosum62.dat";
		
		String cmdString = "./LAKernel/LAkernel_direct " + seq1 + " " + seq2 + " ./LAKernel/" + matrix;
		double res = 0.0;
		
		try {
			
			Process proc = Runtime.getRuntime().exec(cmdString);
		
			proc.getOutputStream().flush();
			proc.getOutputStream().close();
		
			proc.waitFor();
		
		/*
		 * 
		 * extract result
		 * 
		 */

			BufferedReader input = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		
			String line = input.readLine();
			
			res = Double.parseDouble(line);
			
			input.close();
			
			proc.destroy();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while calculating LAkernel.");
		}
		catch (InterruptedException e) {
			System.out.println(e.getMessage());
			System.out.println("InterruptedException while calculating LAkernel.");
		}
		
		return res;
		
	}
	
	public double getNormalizedLAKernel(String seq1, String seq2, String matrix) {
		

		return getLAKernelScore(seq1, seq2, matrix) / Math.sqrt(getLAKernelScore(seq1, seq1, matrix) * getLAKernelScore(seq2, seq2, matrix));
		
		
	}
	
	
	public double getNormalizedMMKernel(String seq1, String seq2, int k, int m, String base_dir) {
		
		
		BufferedWriter bw = null;
		
		String mmkernel_dir = "mismatchkernel/";
		String cmdString = base_dir + mmkernel_dir + "launchMMKernel";
		
		double res = 0.0;
		
		try {
			
		// write input file for kernel
			
			bw = new BufferedWriter(new FileWriter(new File(base_dir + mmkernel_dir + "input.txt")));
			//bw = new BufferedWriter(new FileWriter(new File("./MismatchKernel/data/input.txt")));
			
			bw.write(">seq1\n");
			
			
		// write first sequence	
			
			int len = seq1.length();
			
			String seq = "";
			
			for(int i=0; i < len / 70 + 1; i++) {
				
				for(int j= i * 70; j < Math.min( (i+1) * 70, len ) ; j++) {
					
					seq += seq1.charAt(j) + " ";
					
				}
				
				if(len % 70 > 0) 
					
					seq = seq.substring(0,seq.length()-1) + "\n";
			}
			
			bw.write(seq);
			
			
			bw.write(">seq2\n");
			
			
		// write second sequence	
			
			len = seq2.length();
			
			seq = "";
			
			for(int i=0; i < len / 70 + 1; i++) {
				
				for(int j= i * 70; j < Math.min( (i+1) * 70, len ) ; j++) {
					
					seq += seq2.charAt(j) + " ";
					
				}
				
				if(len % 70 > 0) 
					
					seq = seq.substring(0,seq.length()-1) + "\n";
			}
			
			bw.write(seq);
			
			bw.flush();
			bw.close();
			
			File exec_file = new File(base_dir + mmkernel_dir + "launchMMKernel");
			
			bw = new BufferedWriter(new FileWriter(exec_file));
			
			String rel_dir = " ../../";
			if (base_dir.startsWith("/") || base_dir.startsWith("~")) {
				rel_dir = " ";
			}
			
			bw.write("cd ./MismatchKernel/src/\n");
			bw.write("./string-kernel -K -L " + k + " -D " + m + rel_dir + base_dir + mmkernel_dir + "input.txt\n");
			
			bw.flush();
			bw.close();
			
			exec_file.setExecutable(true);
			
			
		// launch kernel
			
			Process proc = Runtime.getRuntime().exec(cmdString);
		
			proc.getOutputStream().flush();
			proc.getOutputStream().close();
		
			proc.waitFor();
		
		/*
		 * 
		 * extract result
		 * 
		 */
	
			BufferedReader input = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		
			input.readLine();
			
			StringTokenizer strtok = new StringTokenizer(input.readLine());
			
			strtok.nextToken();
			
			double res11 = Double.parseDouble(strtok.nextToken());
			double res12 = Double.parseDouble(strtok.nextToken());
			
			strtok = new StringTokenizer(input.readLine());
			
			strtok.nextToken();
			double res21 = Double.parseDouble(strtok.nextToken());
			double res22 = Double.parseDouble(strtok.nextToken());
			
			if(res12 != res21) {
				System.out.println("Error occurred while calculating kernel value. k(x,y) is not symmetric!");
				System.out.println("k(x,y) = " + res12);
				System.out.println("k(y,x) = " + res21);
				System.out.println("Aborting.");
				System.exit(0);
			}
			
			res = res12 / Math.sqrt(res11 * res22);
			
			input.close();
			
			proc.destroy();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while calculating MMKernel.");
		}
		catch (InterruptedException e) {
			System.out.println(e.getMessage());
			System.out.println("InterruptedException while calculating MMKernel.");
		}
		
		return res;
		
	}

	public static void main(String[] args) throws IOException, BioException {
		
		SequenceAligner aligner = new SequenceAligner("BLOSUM_62.dat", "SW");
		String seq1 = "CCCCCCCCMCGYTGSPEIPQCAGCNQHIVDRFILKVLDRWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW";
		String seq2 = "WWWWWWWWWWWWWWWWMCGYTGSPEPQCACNPHDRFILKVLDRWWWWWWWWWWWWWWWWW";
		String[] res = aligner.getMatchingRegionAndScore(seq1, seq2);
		System.out.println("Alignment score: " + res[0]);
		System.out.println("Start Position:  " + res[1]);
		System.out.println("End Position:    " + res[2]);
	}
  
}

