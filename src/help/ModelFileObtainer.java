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

import java.io.File;

import main.FBPPredictor;
import optimization.MatlignOptimizer;
import optimization.MoStaOptimizer;
import optimization.Optimizer;

public class ModelFileObtainer {
	
	public String obtainModelFile(String class_id) {
		
		return obtainModelFile(FBPPredictor.defaultModelDir, class_id);
	}
	
	public String obtainModelFile(String model_dir, String class_id) {
		
		if (!model_dir.endsWith(File.separator)) {
			model_dir += File.separator;
		}
		
		FileFilter filter = new FileFilter();
		
		filter.setFormat(".*.model");
		
		filter.setDirectory(model_dir + class_id);
		
		String[] files = null;
		if (model_dir.startsWith(File.separator)) {
			files = filter.listFilesFullPath();
		} else {
			files = filter.listFiles();
		}
		
		return model_dir + class_id + "/" + files[0];
	}
}

