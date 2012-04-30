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

package plot;

import java.awt.Color;
import java.awt.Font;

import javax.swing.*;
import org.math.plot.*;
import org.math.plot.plotObjects.BaseLabel;


public class PlotGenerator {
	
	public void plotLine(double[] thresholds, double[] mosta_scores, double[] pred_rates) {
 
		Plot2DPanel plot = new Plot2DPanel();
		plot.addLegend("SOUTH");
		
		BaseLabel plot_title = new BaseLabel("Mean MoSta Scores and Prediction Rates depending on the Threshold", Color.BLUE, 0.5, 1.1);
		plot_title.setFont(new Font("Arial", Font.BOLD, 14));
		plot.addPlotable(plot_title);
		
		// customize axes
		plot.setAxisLabels(new String[] {"Best Match Threshold", "Mean MoSta score / prediction rate"});
		plot.getAxis(0).setLabelPosition(0.5, -0.15);
		plot.getAxis(1).setLabelPosition(-0.15, 0.5);
		plot.getAxis(1).setLabelAngle(-Math.PI / 2);
		
		plot.addLinePlot("Mean MoSta score", thresholds, mosta_scores);
		plot.addLinePlot("Prediction rate", thresholds, pred_rates);

		JFrame frame = new JFrame("SABINE: Prediction rate and accuracy");
		frame.setSize(600, 600);
		frame.setContentPane(plot);
		frame.setVisible(true);

	}
	
	public void plotHist(double[] data, int num_bins, String title) {
 
		Plot2DPanel plot = new Plot2DPanel();
		plot.addHistogramPlot(title, data, num_bins);
		
		// add title
		BaseLabel plot_title = new BaseLabel(title, Color.BLUE, 0.5, 1.1);
		plot_title.setFont(new Font("Arial", Font.BOLD, 14));
		plot.addPlotable(plot_title);

		// customize axes
		plot.setAxisLabels(new String[] {"MoSta score", "Absolute frequency"});
		plot.getAxis(0).setLabelPosition(0.5, -0.15);
		plot.getAxis(1).setLabelPosition(-0.15, 0.5);
		plot.getAxis(1).setLabelAngle(-Math.PI / 2);
		
		// put the PlotPanel in a JFrame
		String frame_title = "SABINE: Distribution of MoSta scores";
		if (title.indexOf("=") != -1) {
			 frame_title += " (t = " + title.substring(title.lastIndexOf("=")+2) + ")";
		}
		JFrame frame = new JFrame(frame_title);
		frame.setSize(600, 600);
		frame.setContentPane(plot);
		frame.setVisible(true);
	}
	
	public void plotBars(double[] data, String title) {
		
		Plot2DPanel plot = new Plot2DPanel();
		plot.addBarPlot(title, data);
		
		//customize axes
		plot.setAxisLabels(new String[] {"model", "MSA / AAE"});
		plot.getAxis(0).setLabelPosition(0.5, -0.15);
		plot.getAxis(1).setLabelPosition(-0.15, 0.5);
		plot.getAxis(1).setLabelAngle(-Math.PI / 2);
		
		JFrame frame = new JFrame("SABINE: My bar plot");
		frame.setSize(600, 600);
		frame.setContentPane(plot);
		frame.setVisible(true);
	}
	
	public static void main(String[] args) {
		//double[] data = new double[] { 0.1,0.1,0.1,0.3,0.1,0.1,0.4,0.4,0.3,0.3,0.2,0.2,0.2,0.2 };
		
		/*
		PlotGenerator plotter = new PlotGenerator();
		plotter.plotLine(new double[] {0.98 , 0.99, 1.0}, 
						 new double[] {0.43, 0.49 , 0.50},
						 new double[] {0.3, 0.23, 0.17});
		*/
		//plotter.plotHist(data, 50, "Best Match Threshold = 0.95");
		
		PlotGenerator plotter = new PlotGenerator();
		plotter.plotBars(new double[] {0.43, 0.743, 0.64}, "test");
	}
	
	
}



