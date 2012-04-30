import javax.annotation.processing.AbstractProcessor;

import de.zbit.util.progressbar.AbstractProgressBar;

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


public class SABINE_Runner implements Runnable {

	private String base_dir;
	private double bmt;
	private boolean dyn_bmt;
	private int mnb;
	private double oft;
	private GUIListener gui;
	private AbstractProgressBar progressBar=null;
	
	public SABINE_Runner(String base_dir, double bmt, boolean dyn_bmt, int mnb, double oft, GUIListener gui) {
		this.base_dir = base_dir;
		this.bmt = bmt;
		this.dyn_bmt = dyn_bmt;
		this.mnb = mnb;
		this.oft = oft;
		this.gui = gui;
	}
  
  public void setProgressBar(AbstractProgressBar progress) {
    this.progressBar = progress;
  }

	public void run() {
		
		// set parameters
		FBPPredictor predictor = new FBPPredictor();
		predictor.dynamic_threshold = dyn_bmt;
		predictor.best_match_threshold = bmt;
		predictor.max_number_of_best_matches = mnb;
		predictor.outlier_filter_threshold = oft;
		predictor.silent = true;
		predictor.gui_output_mode = true;
		predictor.setProgressBar(progressBar);
		

		// create directories for temporary files
		SABINE_Caller dir_creator = new SABINE_Caller();
		//HACK
		dir_creator.createTempDirectories(base_dir);

		// run SABINE on generated input file
		//HACK
		predictor.predictFBP(base_dir + "infile.tmp", base_dir, "trainingsets_public/", null);
		gui.done();
	}

}

