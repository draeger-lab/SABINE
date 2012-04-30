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

public class InfileSplitter {

	
	public void splitInfile(String infile, int num_entries) {
		
		try{
			
			String line;
			int entry_cnt = 0;
			int file_idx = 1;
			int tf_counter = 0;
			int na_counter = 0;
			int line_cnt = 0;
			
			String infile_prefix = infile.substring(0,infile.indexOf('.'));
			String infile_suffix = infile.substring(infile.indexOf('.'));
			
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(infile_prefix + "_" + file_idx + infile_suffix)));

			while ((line = br.readLine()) != null) {
				
				line_cnt++;
				bw.write(line + "\n");
				
				if (entry_cnt == num_entries) {
					
					entry_cnt = 0;
					
					bw.flush();
					bw.close();
					
					bw = new BufferedWriter(new FileWriter(new File(infile_prefix + "_" + ++file_idx + infile_suffix)));
				}
				
				if (line.startsWith("NA  ")) {
					na_counter++;
				}
				
				if (line.startsWith("//")) {
					entry_cnt++;
					tf_counter++;
					
					if (tf_counter != na_counter) {
						System.out.println("Line: " + line_cnt);
						System.exit(0);
					}
				}
			}
			br.close();
			bw.flush();
			bw.close();
		
			System.out.println(tf_counter + " factors processed.");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while calling SABINE.");
		}
	}
	
	
	public static void main(String[] args) {
		
		String infile = "/home/eichner/testdir/missing_mm.input";
		
		InfileSplitter splitter = new InfileSplitter();
		splitter.splitInfile(infile, 10);
	}
}

