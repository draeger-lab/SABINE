package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;


public class ManualAnnotation {

	public HashMap<String, String> annote(File in, HashMap<String, String> hmTransfac) {
		/* 
		 * ?? -- AT_hook
		 * ?? -- Zf_C2H2
		 * ?? -- APSES-type
		 * ?? -- NDT80_PhoG
		 * ?ATF? -- AFT
		 * 
		 * 2.4 Class: Cys6 -- Zn2Cys6
		 * 3.3 Class: Fork head -- Fork-head
		 * 2.2.1 Family: GATA-Factors -- Zf_GATA
		 * 4.4 Class: MADS box -- MADS
		 * 3.4.1 Family: HSF -- HSF_DNA-bind
		 * 4.7 Class: HMG -- HMG_box
		 * 
		*/
		
		int count = 0;
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(in));
			String line = br.readLine();
		
			 while (line != null) {
				 
				 String key = null;
				 String val = null;
				 
				 line = line.trim();
				 
				 if (line.isEmpty() != true) {
					 
					 StringTokenizer st = new StringTokenizer(line);
						
					 val = st.nextToken();
					 key = st.nextToken();
					 
					 //System.out.println(key+" "+val);
				     hmTransfac.put(key.trim(), val.trim());
				     
				     count++;
					 
				 }
				 
				 

				 line = br.readLine();
			 }
			
		 
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		
		System.out.println(count+" manual annotations added to transfac-mapping.");
		return hmTransfac;
		
	}

}
