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
import java.util.StringTokenizer;

public class TFNamePairsFileGenerator {
	
	public void generateTFNamePairsFile(String infile, String outfile) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(infile)));
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			String line = null;
			
			StringTokenizer strtok = null;
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				bw.write(strtok.nextToken() + " ");
				bw.write(strtok.nextToken() + " ");
				bw.write(strtok.nextToken() + "\n");
				
			}
			
			br.close();
			bw.flush();
			bw.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while generating file with TF-namepairs.");
		}
		
	}
	
}

