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

package model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;


public class FeatureSelector {
	
	
	/*
	 * 
	 * selects a precomputed set of features from a given dataset "infile"
	 * 
	 * deletes all other features within this dataset
	 * 
	 * INPUT : ORDERED (!) set of feature indices
	 * 
	 */
	public void selectFeatures(int[] indices, String infile, String outfile) {
		
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(infile)));
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			String line = null;
			
			String feature = null;
			
			StringTokenizer strtok = null;
			
			int index_counter = 0;
			
			
		// go over all data points	
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				index_counter = 0;
				
				
				bw.write(strtok.nextToken());									// label
				
				while(strtok.hasMoreTokens()) {
					
					feature = strtok.nextToken();								// feature
					
					if(feature.startsWith(indices[index_counter] + ":")) {
						
						feature = feature.split(":")[1];
						
						bw.write(" " + (index_counter + 1) + ":" + feature);
						
						index_counter++;
						
						if(index_counter >= indices.length) break;
					}
					
				}
			
				bw.write("\n");
				
			}
			
			bw.flush();
			bw.close();
			br.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while selecting features.");
		}
		
	}

	
	public static void main(String[] args) {
		
		FeatureSelector selector = new FeatureSelector();
		
		selector.selectFeatures(new int[] {1,5,7,9,23}, "data/mosta.lp.24features.att", "data/mosta.lp.5features.att");
		
	}
	
}

