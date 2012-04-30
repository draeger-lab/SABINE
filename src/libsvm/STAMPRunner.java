package libsvm;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;



public class STAMPRunner {
	
	
	/*
	 * please specify the path to your STAMP executable
	 */
	
	String path2STAMP = "./STAMP/";
	public String basedir = "";
	
	public ArrayList<String> runSTAMP(ArrayList<ArrayList<String>> fbps) {
		
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		ArrayList<String> res = new ArrayList<String>();
		
		if (basedir.isEmpty()) {
			basedir = path2STAMP;
		}
		
		try {
				
			
			bw = new BufferedWriter(new FileWriter(new File(basedir + "stamp/input.motifs")));
			
			
		// go over all single fbps	
			
			for(int i=0; i<fbps.size(); i++) {
				
				bw.write("DE\t" + ((10001 + i) + "").substring(1) + "\n");
				
				
			// go over current single fbp
				
				for(int j=0; j<fbps.get(i).size(); j++) {
					
					bw.write(fbps.get(i).get(j) + "\n");
					
				}
				
				bw.write("XX\n");
				
			}	
			
			bw.flush();
			bw.close();
				
			
		// generate shell script

			
			File exec_file = new File(basedir + "stamp/launchSTAMP");
			bw = new BufferedWriter(new FileWriter(exec_file));
			bw.write("export LD_LIBRARY_PATH=" + System.getProperty("user.dir") + "/GSL/gsl-1.10/lib\n");
			bw.write("STAMP/code/STAMP -tf " + basedir + "stamp/input.motifs -sd STAMP/ScoreDists/JaspRand_PCC_SWU.scores -cc PCC -align SWU -ma IR -printpairwise -match STAMP/jaspar.motifs -out " + basedir + "stamp/output\n");

			bw.flush();
			bw.close();
			exec_file.setExecutable(true);
			
			
			

			
		// run STAMP 
					
			run();
	
			
		// extract FBP and return it as ArrayList	
						
			int a,c,g,t,sum;
			
			String max = null;
			
			StringTokenizer strtok = null;
				
			br  = new BufferedReader(new FileReader(new File(basedir + "stamp/outputFBP.txt")));
						
			String line = br.readLine();  					// first line is irrelevant
			
			String fbpline = null;
			
					
			while ((line = br.readLine()) != null) {
						
						
				if(line.startsWith("XX")) {
								
					br.close();
					
					break;
							
				}	
				
				
				
				strtok = new StringTokenizer(line);		
							
				fbpline = strtok.nextToken() + "\t"; // line number
							
							
				a = Math.round(Float.parseFloat(strtok.nextToken()) * 100);
							
				max = "a";
								
								
				c = Math.round(Float.parseFloat(strtok.nextToken()) * 100);
								
				if (c > a) max = "c";
								
								
				g = Math.round(Float.parseFloat(strtok.nextToken()) * 100);
								
				if(g > Math.max(a, c)) max = "g";
								
								
				t = Math.round(Float.parseFloat(strtok.nextToken()) * 100);
							
				if(t > Math.max(g, Math.max(a, c))) max = "t";
								
								
				sum = a + c + g + t;
								
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
							
				if(a + c + g + t != 100) {
					
					System.out.println("Error. Summe: " + ( a + c + g + t ));
				
				}				
				
				fbpline += a + "\t";
				fbpline += c + "\t";
				fbpline += g + "\t";
				fbpline += t + "\t";
							
				fbpline += strtok.nextToken(); // consensus letter
				
				res.add(fbpline);
			
			}
			
			br.close();
			
			
			new File(path2STAMP + "outputFBP.txt").delete();
			
			
		}
		catch(FileNotFoundException fnfe) {
			System.out.println(fnfe.getMessage());
			System.out.println("File not found.");
		}
		catch(IOException ioe) {
			System.out.println("IOException occurred while writing matrices to file.");
		}
		
		return res;
		
	}
	
	
	
	public void run() {
		
		String line = null;
		
		// execute method	
		
		try {
			
			Process proc = Runtime.getRuntime().exec(basedir + "stamp/launchSTAMP");
			
			BufferedReader input = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
		     
			while ((line = input.readLine()) != null) {
				
		        System.out.println(line);
		    
			}
			
			input.close();
			
			proc.waitFor();
			
			
			
		}
		catch(IOException e) {
			System.out.println("IOException while executing Script!");
			e.printStackTrace();
		} 
		catch (InterruptedException e) {
			System.out.println("InterruptedException while executing Perl Script!");
		}
		
	}
	
}
