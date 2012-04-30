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

import java.io.BufferedWriter;
import java.util.ArrayList;

public interface SABINEParser {
	
	public ArrayList<String> get_tf_names();
	public ArrayList<String> get_sequences1();
	public ArrayList<String> get_sequences2();
	public ArrayList<String> get_species();
	public ArrayList<ArrayList<String>> get_domains();
	public ArrayList<ArrayList<String>> get_matrices();
	
	public void parseAll(int class_id, String inputfile, String labelfile);
	public void writeInputFile(BufferedWriter bw, int[] relevant_factors, int class_id);
}

