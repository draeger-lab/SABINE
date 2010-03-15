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

import help.FileFilter;

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
 * eliminates numerical instabilities, which are sometimes caused by 3
 *  
 * of the 12 used substitution matrices
 * 
 */
public class FeatureRepairer {
	
	public boolean silent = false;
	
	
	public void repairAllFeatureFiles(String indir) {
		
		FileFilter filter = new FileFilter();
		
		filter.setFormat(".*.out");
		
		filter.setDirectory(indir);
		
		String[] files;
		if (indir.startsWith("/")) {
			files = filter.listFilesFullPath();
		}
		else {
			files = filter.listFiles();
		}
		
		if (! silent) System.out.println();
		
		for(int i=0; i<files.length; i++) {
			
			if (! silent) System.out.print("  Screening " + files[i] + "...");
			
			repairFeatureFile(indir + "/" + files[i], indir + "/" + files[i]);
			
			if (! silent) System.out.println("done.");
			
		}
		
		
	}
	
	
	public void repairFeatureFile(String infile, String outfile) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		 
		ArrayList<String> entries = new ArrayList<String>();
	
		try {
			 
			 br = new BufferedReader(new FileReader(new File(infile)));
			 
			 String line = null;
			 
			 StringTokenizer strtok = null;
			 
			 
			 boolean repaired = false;
			
			 
			 while((line = br.readLine()) != null) {
				 
				 strtok = new StringTokenizer(line);
				 
				 String entry = "";
				 
				 entry += strtok.nextToken() + " ";
				 
				 entry += strtok.nextToken() + " ";
				 
				 entry += strtok.nextToken() + " ";
				 
				 entry += strtok.nextToken() + " ";
				 
				 double score = Double.parseDouble(strtok.nextToken());
				 
				 if(score == Double.POSITIVE_INFINITY) { score =  100.0; repaired = true; }
				 if(score == Double.NEGATIVE_INFINITY) { score = -100.0; repaired = true; }
				 if(score == Double.NaN) 			   { score = 	0.0; repaired = true; }
				 
				 
				 entries.add(entry + score);
				 
			 }
			 
			 if(repaired && !silent) System.out.print("repairing...");
			 
			 br.close();
			 
			 
			 bw = new BufferedWriter(new FileWriter(new File(outfile)));
			 
			 
			 for(int i=0; i<entries.size(); i++) {
				 
				 bw.write(entries.get(i) + "\n");
				 
			 }
			 
			 bw.flush();
			 bw.close();
			 
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while repairing features.");
		}
		
	}
	
	public static void main(String[] args) {
		
		String test_dir = "/rahome/eichner/sabine_project/SABINE_instances/basedir1/relevantpairs";
		//String test_dir = "basedir1/relevantpairs";
		
		FeatureRepairer repairer = new FeatureRepairer();
		repairer.repairAllFeatureFiles(test_dir);
	}
	
}

