/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */
package util.parsers.tools.IPRscan;
/*
 * Converter Class
 */
import java.util.ArrayList;

public class Converter {
	
	
	// convert ArrayList to Array of booleans
	public boolean[] l2r(ArrayList<Boolean> flags_arr) {
		boolean[] flags = new boolean[flags_arr.size()];
		
		for (int i = 0; i < flags_arr.size(); i++) {
			if (flags_arr.get(i)) {
				flags[i] = true;
			}
		}
		return flags;
	}
	
	// convert absolute MA values to relative MA values
	public ArrayList<ArrayList<ArrayList<String>>> a2r(ArrayList<ArrayList<ArrayList<String>>> matrices_in) {
		ArrayList<ArrayList<ArrayList<String>>> matrices_out = new ArrayList<ArrayList<ArrayList<String>>>();
		ArrayList<ArrayList<String>> curr_pfms = null;
		ArrayList<String> curr_pfm = null;
		
		
		for (int i = 0; i < matrices_in.size(); i++) {
			curr_pfms = new ArrayList<ArrayList<String>>();
			for (int j = 0; j < matrices_in.get(i).size(); j++) {
				curr_pfm = new ArrayList<String>();
				
				
				String[] aarr = matrices_in.get(i).get(j).get(0).split("   ");
				String[] carr = matrices_in.get(i).get(j).get(1).split("   ");
				String[] garr = matrices_in.get(i).get(j).get(2).split("   ");
				String[] tarr = matrices_in.get(i).get(j).get(3).split("   ");
				
				String astr = "";
				String cstr = "";
				String gstr = "";
				String tstr = "";
				
				for (int k = 0; k < aarr.length; k++) {
					
					float a = Float.valueOf(aarr[k]).floatValue();
					float c = Float.valueOf(carr[k]).floatValue();
					float g = Float.valueOf(garr[k]).floatValue();
					float t = Float.valueOf(tarr[k]).floatValue();
					
					float sum = a+c+g+t;
					
					a/=sum;
					c/=sum;
					g/=sum;
					t/=sum;

					astr = astr.concat(a + "   ");
					cstr = cstr.concat(c + "   ");
					gstr = gstr.concat(g + "   ");
					tstr = tstr.concat(t + "   ");
				}
				
				curr_pfm.add(astr);
				curr_pfm.add(cstr);
				curr_pfm.add(gstr);
				curr_pfm.add(tstr);
				
				curr_pfms.add(curr_pfm);
			}
			matrices_out.add(curr_pfms);
		}
		return matrices_out;
	}

}
