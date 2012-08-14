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


import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.StringTokenizer;

import libsvm.PWMOutlierFilter;
import libsvm.STAMPRunner;

import optimization.MatrixConverter;
import optimization.MoStaOptimizer;



public class PFMFormatConverter {
	
	public String basedir = "";
	public boolean silent = false;
	
	
	public ArrayList<ArrayList<String>> convertAllTransfacToSTAMP(ArrayList<String[]> transfac_pfms) {
		
		ArrayList<String> curr_pfm;
		ArrayList<ArrayList<String>> all_pfms = new ArrayList<ArrayList<String>>();
		
		for (int i=0; i<transfac_pfms.size(); i++) {
			curr_pfm = convertTransfacToSTAMP(transfac_pfms.get(i));
			all_pfms.add(curr_pfm);
		}
		return all_pfms;
	}
	
	/*
	 *  converts TRANSFAC format to STAMP format
	 */
	
	public ArrayList<String> convertTransfacToSTAMP(String[] transfac_pfm) {
		
		ArrayList<String> stamp_pfm = new ArrayList<String>();
		String line, max, cons_letter;
		float a, c, g, t;
		int sum;
		MatrixConverter converter = new MatrixConverter(); 
		
		StringTokenizer strtok1 = new StringTokenizer(transfac_pfm[0]);
		StringTokenizer strtok2 = new StringTokenizer(transfac_pfm[1]);
		StringTokenizer strtok3 = new StringTokenizer(transfac_pfm[2]);
		StringTokenizer strtok4 = new StringTokenizer(transfac_pfm[3]);
		
		int num_tokens = strtok1.countTokens();
		
		for (int i=0; i<num_tokens; i++) {
			
			a = Float.parseFloat(strtok1.nextToken().trim());
			c = Float.parseFloat(strtok2.nextToken().trim());
			g = Float.parseFloat(strtok3.nextToken().trim());
			t = Float.parseFloat(strtok4.nextToken().trim());
			
			cons_letter = converter.calculateConsensusLetter(a, c, g, t);
			
			a = Math.round(a * 100);
			c = Math.round(c * 100);
			g = Math.round(g * 100);
			t = Math.round(t * 100);
			
			max = "a";				
			if (c > a) max = "c";		
			if(g > Math.max(a, c)) max = "g";			
			if(t > Math.max(g, Math.max(a, c))) max = "t";
									
			sum = (int) (a + c + g + t);
							
			if(sum > 100) {				
				for(int j=100; j<sum; j++) {
									
					if(max.equals("a")) a--;
					if(max.equals("c")) c--;
					if(max.equals("g")) g--;
					if(max.equals("t")) t--;					
				}			
			}			
			else {				
				if(sum < 100) {				
					for(int j=sum; j<100; j++) {
									
						if(max.equals("a")) a++;
						if(max.equals("c")) c++;
						if(max.equals("g")) g++;
						if(max.equals("t")) t++;			
					}				
				}
			}
			line = i + "\t" + (int) a + "\t" + (int) c + "\t" + (int) g + "\t" + (int) t + "\t" + cons_letter;
			stamp_pfm.add(line);	
		}
		return stamp_pfm;
	}
	
	public String[] convertSTAMPToTransfac(ArrayList<String> stamp_pfm) {
	
		StringTokenizer strtok = new StringTokenizer(stamp_pfm.get(0));;
		
		DecimalFormat fmt = new DecimalFormat();
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
        symbs.setDecimalSeparator('.');
        fmt.setDecimalFormatSymbols(symbs);
        fmt.setMaximumFractionDigits(4);
        fmt.setMinimumFractionDigits(4);
		
		strtok.nextToken();
		float a = Float.parseFloat(strtok.nextToken().trim());
		float c = Float.parseFloat(strtok.nextToken().trim());
		float g	= Float.parseFloat(strtok.nextToken().trim());
		float t	= Float.parseFloat(strtok.nextToken().trim());
		float sum = a + c + g + t;
		
		String line1 = "" + fmt.format(a / sum);
		String line2 = "" + fmt.format(c / sum);
		String line3 = "" + fmt.format(g / sum);
		String line4 = "" + fmt.format(t / sum);
		
		for (int i=1; i<stamp_pfm.size(); i++) {
			strtok = new StringTokenizer(stamp_pfm.get(i));
			
			strtok.nextToken();                    			  // Index
			a = Float.parseFloat(strtok.nextToken().trim());  // A
			c = Float.parseFloat(strtok.nextToken().trim());  // C
			g = Float.parseFloat(strtok.nextToken().trim());  // G
			t = Float.parseFloat(strtok.nextToken().trim());  // T
			sum = a + c + g + t;
			
			line1 += "   " + fmt.format(a / sum);
			line2 += "   " + fmt.format(c / sum);
			line3 += "   " + fmt.format(g / sum);
			line4 += "   " + fmt.format(t / sum);
			
		}
	    return new String[] {line1, line2, line3, line4};
	}
	
	public  ArrayList<ArrayList<String>> mergePFMnames(ArrayList<ArrayList<String>> pwm_names) {
		
		ArrayList<ArrayList<String>> merged_pwm_names = new ArrayList<ArrayList<String>>();
		ArrayList<String> merged_pwm_name = new ArrayList<String>();
		String curr_pwm_name;
		
		for (int i=0; i<pwm_names.size(); i++) {
			
			merged_pwm_name = new ArrayList<String>();
			curr_pwm_name = pwm_names.get(i).get(0);
			
			for (int j=1; j<pwm_names.get(i).size(); j++) {
				curr_pwm_name += "," + pwm_names.get(i).get(j);
			}
			merged_pwm_name.add(curr_pwm_name);
			merged_pwm_names.add(merged_pwm_name);
		}
		return merged_pwm_names;
	}
	
	public ArrayList<String> mergePFMs(ArrayList<ArrayList<String>> pwms) {
		
		ArrayList<String> fbp = new ArrayList<String>(); 
		
		// filter outlier PWMs according to MoSta score
		PWMOutlierFilter filter = new PWMOutlierFilter();
		filter.silent = silent;
		
		MoStaOptimizer  opt = new MoStaOptimizer();
		opt.basedir = basedir;
		pwms = filter.filterOutliers(pwms, opt);
		
		// merge remaining PWMs using STAMP
		STAMPRunner stamp_runner = new STAMPRunner();
		stamp_runner.basedir = basedir;
		fbp = stamp_runner.runSTAMP(pwms);
		
		return fbp;
	}
	
	public static void main(String[] args) {
		
	}

}

