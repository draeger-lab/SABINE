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
import java.io.File;
import java.io.FileNotFoundException;
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
	 * fields
	 */
	
	
	FiniteAlphabet alphabet = null;
    
	SubstitutionMatrix matrix = null;
	
	SequenceAlignment alignment = null;
	
	StringTokenizer strtok = null;
	
	/*
	 * constructors
	 */
	
	public SequenceAligner(String matrixfile, String type) throws IOException, BioException {
		
		alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
		
		
		matrix = new SubstitutionMatrix(alphabet, new File("substitutionMatrices" + File.separator + matrixfile));
		
		if(type.equals("NW")) {
			alignment = new NeedlemanWunsch( 
			        0, 		// match
			        0,		// replace
			        11,  	// insert
			        11,		// delete
			        1,  	// gapExtend
			        matrix 	// SubstitutionMatrix
			      );
			      
		}
		
		if(type.equals("NW_UNGAPPED")) {
			alignment = new NeedlemanWunsch( 
			        0, 						// match
			        3,						// replace
			        Integer.MAX_VALUE,  	// insert
			        Integer.MAX_VALUE,		// delete
			        Integer.MAX_VALUE,  	// gapExtend
			        matrix 					// SubstitutionMatrix
			      );
			      
		}
		
		if(type.equals("SW")) {
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
	
	
	
	public SequenceAligner(int match, int mismatch, String type) throws IOException, BioException {
		
		alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
		
		matrix = new SubstitutionMatrix(alphabet, match, mismatch);
		
		if(type.equals("NW")) {
			alignment = new NeedlemanWunsch( 
			        -1, 	// match
			        1,		// replace
			        1,  	// insert
			        1,		// delete
			        1,  	// gapExtend
			        matrix 					// SubstitutionMatrix
			      );
			      
		}
		else {
			alignment = new SmithWaterman(
			        0,      // match
			        1,      // replace 
			        1,      // insert
			        1,      // delete
			        1,      // gapExtend
			        matrix  				// SubstitutionMatrix
			      );
		}
	
	}
	
	/*
	 * methods
	 */
	
	
	public double getSequenceSimilarity(String seq1, String seq2) {
	
		double similarity = 0.0;
		
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
		
		String upperseq = null;
		
	// goto start of alignment	
		for(int i=0; i<6; i++) line = strtok1.nextToken();
		
	//	System.out.println("Line:" + line);
		
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
			
			line = line.substring(line.length()- upperseq.length() - 1,line.length()-1); // trim whitespaces before |'s
			line = line.substring(0,upperseq.length()); // trim whitespace after |'s
			
	//		System.out.println("Line:  " + line + " , Length: " + line.length());
			
			
			for(int i=0; i<line.length(); i++) {
				
				count_all++;
				
				if(line.charAt(i) == '|') 
					count_match++;
			}
			
		if(strtok1.hasMoreTokens())	line = strtok1.nextToken(); // goto next line
		if(strtok1.hasMoreTokens())	line = strtok1.nextToken(); // goto next line
		if(strtok1.hasMoreTokens())	line = strtok1.nextToken(); // goto next line
		
			
		}
		
		similarity = (double) count_match / (double) count_all;
		
	//	System.out.println("\n#Matches:   " + count_match);
	//	System.out.println("#Positions: " + count_all);
	//	System.out.println("Similarity: " + fmt.format(similarity) + "%");
		
		return similarity;
	}
	
	public double getAlignmentScore(String seq1, String seq2) throws Exception {
	      
		  Sequence query  = ProteinTools.createProteinSequence(seq1, "query");
	      Sequence target = ProteinTools.createProteinSequence(seq2, "target");
	      
	      alignment.pairwiseAlignment(query, target);
	      
	      strtok = new StringTokenizer(alignment.getAlignmentString());
	      
	      while(!strtok.nextToken().equals("Score:"));
	      
	   //   System.out.println(alignment.getAlignmentString());
	      
	      
		  return Double.parseDouble(strtok.nextToken());
	      
	}
	
	/*
	 * 
	 * calculates a sequence identity which ignores mismatches whose SM-Scores are above a given threshold
	 * 
	 */
	
	public static double getSMBasedIdentity(String seq1, String seq2, String matrix, double threshold) {
	
		SequenceAligner aligner = null;
		
		double similarity = 0.0;
		
		try {
			
			aligner = new SequenceAligner(matrix, "NW");
		
		}
		catch (BioException be) {
			System.out.println(be.getMessage());
		}
		catch(FileNotFoundException fnfe) {
			System.out.println(fnfe.getMessage());
		}
		catch (IOException ioe) {
			System.out.println(ioe.getMessage());
		}
		
		
		try {
			
			/*
			 * 
			 * align sequences and calculate score
			 * 
			 */
			
			
			
			Sequence query  = ProteinTools.createProteinSequence(seq1, "query");
		    Sequence target = ProteinTools.createProteinSequence(seq2, "target");
		      
		    aligner.alignment.pairwiseAlignment(query, target);
		
		    
		    
		//    System.out.println(aligner.alignment.getAlignmentString());
		    
		    
		    
		    StringTokenizer strtok  = new StringTokenizer(aligner.alignment.getAlignmentString(), "\n");
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
				
				
	//			System.out.println(upperseq + "\n" + dashes + "\n" + lowerseq + "\n");
				
				
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
							
						aligner.matrix.getValueAt(upperSymbol, lowerSymbol) >= threshold ) ) 
						
						
						count_match++;
				}
				
			
				if(strtok.hasMoreTokens())	line = strtok.nextToken(); // goto next line
				if(strtok.hasMoreTokens())	line = strtok.nextToken(); // goto next line
			
			} 
			
			similarity = (double) count_match / (double) count_all;
		
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
	 * 
	 * calculates the substitution matrix based sequence similarity of two sequences
	 * 
	 */
	
	public static double getSMBasedSimilarity(String seq1, String seq2, String matrix) {
		
		SequenceAligner aligner = null;
		
		double similarity = 0.0;
		
		double self_score1 = 0.0;
		double self_score2 = 0.0;
		
		StringTokenizer strtok1 = null;
		
		try {
		
			aligner = new SequenceAligner(matrix, "NW");
		
		}
		catch (BioException be) {
			System.out.println(be.getMessage());
		}
		catch(FileNotFoundException fnfe) {
			System.out.println(fnfe.getMessage());
		}
		catch (IOException ioe) {
			System.out.println(ioe.getMessage());
		}
		
		
		
		try {
			
			/*
			 * 
			 * align sequences and calculate score
			 * 
			 */
			
			Sequence query  = ProteinTools.createProteinSequence(seq1, "query");
		    Sequence target = ProteinTools.createProteinSequence(seq2, "target");
		      
		    aligner.alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(aligner.alignment.getAlignmentString());
			
			
			
		//	System.out.println(aligner.alignment.getAlignmentString());
			
			
			
			
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
		      
		    aligner.alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(aligner.alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			self_score1 = Double.parseDouble(strtok1.nextToken());
			
			
		// compare seq2 vs. itself	
			
			query  = ProteinTools.createProteinSequence(seq2, "query");
		    target = ProteinTools.createProteinSequence(seq2, "target");
		      
		    aligner.alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(aligner.alignment.getAlignmentString());
			
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
	
	
	
	/*
	 * 
	 * calculates the substitution matrix based sequence similarity of two sequences
	 * 
	 */
	
	public static double getSMBasedSecondaryStructureSimilarity(String seq1, String seq2, String matrix) {
		
		SequenceAligner aligner = null;
		
		double similarity = 0.0;
		
		double self_score1 = 0.0;
		double self_score2 = 0.0;
		
		StringTokenizer strtok1 = null;
		
		try {
		
			aligner = new SequenceAligner("DNA", matrix, "NW");
		
		}
		catch (BioException be) {
			System.out.println(be.getMessage());
		}
		catch(FileNotFoundException fnfe) {
			System.out.println(fnfe.getMessage());
		}
		catch (IOException ioe) {
			System.out.println(ioe.getMessage());
		}
		
		
		
		try {
			
			/*
			 * 
			 * align sequences and calculate score
			 * 
			 */
			
			Sequence query  = DNATools.createDNASequence(seq1, "query");
		    Sequence target = DNATools.createDNASequence(seq2, "target");
		      
		    aligner.alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(aligner.alignment.getAlignmentString());
			
			
			
	//		System.out.println(aligner.alignment.getAlignmentString());
			
			
			
			
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
		      
		    aligner.alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(aligner.alignment.getAlignmentString());
			
			while(!strtok1.nextToken().equals("Score:"));
			
			self_score1 = Double.parseDouble(strtok1.nextToken());
			
			
		// compare seq2 vs. itself	
			
			query  = DNATools.createDNASequence(seq2, "query");
		    target = DNATools.createDNASequence(seq2, "target");
		      
		    aligner.alignment.pairwiseAlignment(query, target);
		
			strtok1 = new StringTokenizer(aligner.alignment.getAlignmentString());
			
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
	
	
	
	/*
	 * inefficient method for getting the SM-Score of two AAs
	 */
	public static double getAlignmentScore(String char1, String char2, String matrix) {
		
		
		double res = 0.0;
		
		if(char1.length() + char2.length() != 2) {
			System.out.println("Function can only handle single characters.");
			System.exit(0);
		}
		
		SequenceAligner aligner = null;
		
		try {
			
			aligner = new SequenceAligner(matrix, "NW_UNGAPPED");
		
		}
		catch (BioException be) {
			System.out.println(be.getMessage());
		}
		catch(FileNotFoundException fnfe) {
			System.out.println(fnfe.getMessage());
		}
		catch (IOException ioe) {
			System.out.println(ioe.getMessage());
		}
		
		try {
			
			Sequence query  = ProteinTools.createProteinSequence(char1, "query");
		    Sequence target = ProteinTools.createProteinSequence(char2, "target");
		      
		    aligner.alignment.pairwiseAlignment(query, target);
		
			StringTokenizer strtok = new StringTokenizer(aligner.alignment.getAlignmentString());
			
			while(!strtok.nextToken().equals("Score:"));
			
			res = Double.parseDouble(strtok.nextToken());
		}
		catch(BioException be) {
			System.out.println(be.getMessage());
		} 
		catch (Exception e) {
			e.getMessage();
			e.printStackTrace();
		}
		
		return res;
		
	}
	
	public static Symbol getProteinSymbol(char c) {
		
		Symbol res = null;
		
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
	
	
	public static double getLAKernelScore(String seq1, String seq2) {
		
		
		String cmdString = "./LAKernel/LAkernel_direct " + seq1 + " " + seq2 + " ./LAKernel/blosum62.dat";
		
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
	
	public static double getNormalizedLAKernel(String seq1, String seq2) {
		

		return getLAKernelScore(seq1, seq2) / Math.sqrt(getLAKernelScore(seq1, seq1) * getLAKernelScore(seq2, seq2));
		
		
	}
	
	
	
	
  public static void main (String args[]) {
    
    try {
    	
    	//YICKECQRKFSSGHHLTRHKKSVH
    	//HSCPKCGKRFKRRDHVLQHLNKKIPC
    	//YICRECNRQFSSGHHLTRHKKSVH
    	//HSCPRCGKRFKRRDHVLQHLNKKIPC
    	
    	String seq1 = "LDDGYRWRKYGQKVVKGNPNPRSYYKCTQVGCPVRKHVERASHDLRAVITTYEGKHNHDV";
    	String seq2 = "LDDGYRWRKYGQKVVKGNPNPRSYYKCTSQGCPVRKHVERASHDIRSVITTYEGKHNHDV";
    	
    	System.out.println("Similarity: " + getSMBasedSimilarity(seq1, seq2, "BLOSUM_62.dat"));
    	
    //	System.out.println("Alignment Score: " + getSMBasedSecondaryStructureSimilarity(seq1, seq2, "SS-matrix.dat"));
    	
    //	SequenceAligner aligner = new SequenceAligner(1, -1, "NW");
    	
    //	System.out.println("\nSimilarity: " + aligner.getSequenceSimilarity(seq1, seq2));
    	
    } catch (Exception e) {
  
    	e.printStackTrace();
    
    }
    
  }
  
}

