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
import java.io.InputStreamReader;

import resources.Resource;

public class FileCopier {
	
	public static void copy(String infile, String outfile) {
		copy(infile, outfile, false);
	}
			
	public static void copy(String infile, String outfile, boolean isResource) {
		
		try {
			BufferedReader br = null;
			if (isResource) {
				br = new BufferedReader(new InputStreamReader(Resource.class.getResourceAsStream(infile)));
			} else {
				br = new BufferedReader(new FileReader(new File(infile)));
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			
			String line;
			while ((line = br.readLine()) != null) {
				bw.write(line + "\n");
			}
			br.close();
			bw.flush();
			bw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		copy("/home/eichner/test.txt", "/home/eichner/test1.txt");
	}

}

