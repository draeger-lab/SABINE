package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;


public class PWMParser {
	
	ArrayList<String> PWMnames = new ArrayList<String>();
	ArrayList<ArrayList<String>> pwms = new ArrayList<ArrayList<String>>();
	
	//MN
	ArrayList<String> MNs = new ArrayList<String>();
	
	
	public void parse(ArrayList<ZipEntry> ZipEntries, File in) {
		
		ZipFile zf = null;
		try {
			zf = new ZipFile(in);

	    Iterator<ZipEntry> itr0 = ZipEntries.iterator();
	    while (itr0.hasNext()) {
	     
	    	ZipEntry element = itr0.next();

	    	BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(element)));
	    	
	    	String line = null;
	    	
	    	// MN
	    	MNs.add(element.toString());
	    		    	
	    	// read first line
	    	line = br.readLine();
	    	// add filtered names
	    	PWMnames.add(filterName(line));
	    		    	
	    	ArrayList<String> tmp = new ArrayList<String>();
	    	
	    	while ((line = br.readLine()) != null) {

	    		if (line.isEmpty()!=true) {
	    			
	    			tmp.add(line.trim());
	    			
	    		}

			}
	    	// add pwm-entries to corresponding name
	    	if (tmp.size() == 4) {
	    		pwms.add(tmp);
	    	} else {
	    		System.out.println("ERROR parsing PWM!");
	    	}
	    	

			br.close();
  
	    }
	    
		} catch (ZipException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		removedubs();

	}
	
	
	// filter name out of line
	private String filterName(String line) {
		String name = null;
		
		StringTokenizer st = new StringTokenizer(line);
		
		if (st.countTokens() == 1) {
			
			Pattern re = Pattern.compile("_\\d");
			
			Matcher ma = re.matcher(line);
						
			if (ma.find()) {
				
				String[] cln_name = line.split("_\\d");
				name = cln_name[0];
			}
		}
		
		String tmp_name = st.nextToken();
		
		if (line.startsWith("Gene:") || line.startsWith("Protein:")) {
			
			tmp_name = st.nextToken();

			Pattern re = Pattern.compile("-\\D");
			
			Matcher ma = re.matcher(tmp_name);
			
			if (ma.find()) {
				
				String[] cln_name = tmp_name.split("-\\D");
				name = cln_name[0];
				
			} else {
				
				name = tmp_name;
				
			}
		
		}
		
		if (name == null) {
			System.out.println("ERROR parsing PWMs-names!");
		}
		
		return name;
	}
	
	ArrayList<String> PWMnames_dub = new ArrayList<String>();
	ArrayList<ArrayList<String>> pwms_dub = new ArrayList<ArrayList<String>>();
	
	// MN
	ArrayList<String> MNs_dub = new ArrayList<String>();
	
	
	// check for dubs
	public void removedubs() {
		
		System.out.println(PWMnames.size()+" PWMs parsed.");
		
		for (int i = 0; i < PWMnames.size(); i++) {
			
				if ((PWMnames_dub.contains(PWMnames.get(i))) ) {
					
					//System.out.println(PWMnames.get(i)+" "+pwms.get(i));
					
				} else {
										
					PWMnames_dub.add(PWMnames.get(i));
					pwms_dub.add(pwms.get(i));
					
					//MN
					MNs_dub.add(MNs.get(i));
				}
	    	  	    	      	  
		}
		
		System.out.println(PWMnames.size()-PWMnames_dub.size()+" duplicate PWMs removed.");

    }
	
	

	
	// getter...
	public ArrayList<String> getPWMnames() {
		return PWMnames_dub;
	}

	public ArrayList<ArrayList<String>> getPwms() {
		return pwms_dub;
	}
	
	public ArrayList<String> getMNs() {
		return MNs_dub;
	}
	

}
