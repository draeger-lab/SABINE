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

import java.util.ArrayList;


import model.LabelFileGenerator;

public class MatrixFileParser {

	ArrayList<String> matrix_names;
	ArrayList<ArrayList<String>> matrices;
	
	public void readMatrices(String infile) {
		
		LabelFileGenerator matrix_parser = new LabelFileGenerator();
		matrix_names = matrix_parser.getFactorNames(infile);
		matrices = matrix_parser.getAllFBPs(infile);
	}
	
	public ArrayList<String> obtainMatrix(String factor_name) {
		int idx = matrix_names.indexOf(factor_name);
		return matrices.get(idx);
	}
	
	
	
	public static void main(String[] args) {
		
		MatrixFileParser matrixfile_parser = new MatrixFileParser();
		matrixfile_parser.readMatrices("/home/jei/SABINE/trainingsets/trainingsset_all_classes.pfms");
	}

}

