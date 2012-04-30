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
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

public class IrrelevantPairIdentifier {
	
	public boolean[] identifyIrrelevantPairs(String infile, double relevance_threshold, int numPairs) {
		
		
		boolean[] irrelevantPairs = new boolean[numPairs];
		
		String line = null;
		
		StringTokenizer strtok = null;
		
		BufferedReader br = null;
		 
		try {
			 
			 br = new BufferedReader(new FileReader(new File(infile)));
			 
			 int line_count = 0;
			 
			 
			 while((line = br.readLine()) != null) {
			
				 
				 strtok = new StringTokenizer(line);
				 
				 
				 strtok.nextToken();				// query
				 
				 strtok.nextToken();				// vs.
				 
				 strtok.nextToken();				// tf of training set
				 
				 strtok.nextToken();				// :
				 
				 
				 if(Double.parseDouble(strtok.nextToken()) < relevance_threshold) {
					 
					 irrelevantPairs[line_count] = true;
					 
				 }
				 
				 else {
					 
					 irrelevantPairs[line_count] = false;
					 
				 }
				 
				 line_count++;
				 
			 }
			 
			 br.close();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while identifying relevant tf pairs.");
		}
		
		
		return irrelevantPairs;
		
		
	}
	
}

