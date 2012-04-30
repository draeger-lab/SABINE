package libsvm;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;

import optimization.Optimizer;


public class PWMOutlierFilter {
	
	
	public double tolerance = 0.5;
	public boolean silent = false;
	
	
	public ArrayList<ArrayList<String>> filterOutliers(ArrayList<ArrayList<String>> matrices, Optimizer optimizer) {
		
		DecimalFormat fmt = new DecimalFormat();
		
		fmt.setMaximumFractionDigits(2);
		fmt.setMinimumFractionDigits(2);
			
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
		symbs.setDecimalSeparator('.');
		fmt.setDecimalFormatSymbols(symbs);
		
		
		ArrayList<ArrayList<String>> res = new ArrayList<ArrayList<String>>();
		
		
		
			
			boolean[] isOutlier = null;
				
			// calculate all pairwise scores for matrices of this TF
				
				if(matrices.size() > 2) {
				
					double[][] distances = new double[matrices.size()][matrices.size()];
					
					for(int i=0; i<matrices.size(); i++) {
						
						for(int j=i; j<matrices.size(); j++) {
							
							distances[i][j] = optimizer.compareMatrices(matrices.get(i), matrices.get(j), null);
							
						}
						
					}
						
					
				// normalize pairwise scores	
					
					for(int i=0; i<matrices.size()-1; i++) {
						
						for(int j=i+1; j<matrices.size(); j++) {
							
							distances[i][j] = distances[i][j] / Math.sqrt(distances[i][i] * distances[j][j]);
							distances[j][i] = distances[i][j];
						}
					
					}
					
					for(int i=0; i<matrices.size(); i++) {
						
						distances[i][i] = 1.0;
						
					}
					
					
				// transform scores into distances
					
					for(int i=0; i<matrices.size(); i++) {
						
						for(int j=i+1; j<matrices.size(); j++) {
							
							distances[i][j] = 1-distances[i][j];
							distances[j][i] = 1-distances[j][i];
						}
						
						distances[i][i] = 1-distances[i][i];
	
					}
					
					
					
				// print distance matrix	
				
					if (! silent) {	
					
						System.out.println("    Distance Matrix:\n");
					
						for(int i=0; i<distances.length; i++) {
						
							System.out.print("    ");
						
							for(int j=0; j<distances[i].length; j++) {
							
								System.out.print(fmt.format(distances[i][j]) + "  ");
							
							}
						
							System.out.println();
						
						}
					
						System.out.println();
					}	
					
					
					
				// calculate mean pairwise distance
					
					double meanDist = 0;
					double numDists  = 0;
					
					for(int i=0; i<matrices.size()-1; i++) {
						
						for(int j=i+1; j<matrices.size(); j++) {
							
							meanDist += distances[i][j];
							numDists++;
							
						}
	
					}
					
					meanDist /= numDists;
					
					
					if(numDists != matrices.size() * (matrices.size() - 1) * 0.5) {
						
						System.out.println("Fatal Error. Aborting.");
						System.exit(0);
						
					}
					
					if (! silent) {
						System.out.println("    Mean Distance:\t\t\t" + fmt.format(meanDist) + "\n");
					}
					
				// add pseudocount in order to maintain numerical stability	
					
					if(meanDist == 0.0) {
						
						meanDist = meanDist + 0.0000000000001;
						
					}
					
					isOutlier = new boolean[matrices.size()];
					
					
				// identify outliers
					
					for(int i=0; i<matrices.size(); i++) {
						
						double currMeanDist = 0.0;
						
						for(int j=0; j<matrices.size(); j++) {
							
							currMeanDist += distances[i][j];
							
						}
						
						currMeanDist /= matrices.size()-1;
						
						if (! silent) {
							System.out.print("    Mean Distance of Matrix " + i + ":\t" + fmt.format(currMeanDist) + " (" + fmt.format(currMeanDist / meanDist) + ")");
						}
							
						if(currMeanDist > (1.0 + tolerance) * meanDist) {
							
							isOutlier[i] = true;
							
							if (! silent) {
								System.out.println(" --> Outlier!");
							}
						}
						
						else {
							
							isOutlier[i] = false;
							
							if (! silent) {
								System.out.println(" --> OK!");
							}
						}
						
					}
					if (! silent) {
						System.out.println();
					}
				} // endif: more than 2 matrices for this TF were found --> filter was applied	
				
				
			// TF entry contains only one or two matrices --> filtering is not necessary	
				
				else {
					
					isOutlier = new boolean[matrices.size()];
					
					for(int i=0; i<isOutlier.length; i++) {
					
						isOutlier[i] = false;
					
					}
					
				}
				
				
			// add non-outlier matrices of this TF to result 
				
				for(int i=0; i<matrices.size(); i++) {
					
					
					if(!isOutlier[i]) {
						
						res.add(matrices.get(i));
						
					}
					
				}
							
				return res;
	}	
	
}
