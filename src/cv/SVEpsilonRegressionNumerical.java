package cv;

import libsvm.svm_train;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;

public class SVEpsilonRegressionNumerical {

	/**
	 * 
	 * epsilon regression
	 * 
	 * @param str_files
	 * @param str_dir
	 */
	
	public static boolean silent = false;
	
	public static void performEpsilonRegression(String str_file, String str_dir, FileWriter fw) {

		// grid search for C
		for (int log2c = CVGlobalSettings.log2c_start; log2c <= CVGlobalSettings.log2c_end; log2c = log2c + CVGlobalSettings.log2c_incr) {
 
			// grid search for epsilon
			for (int log2eps = CVGlobalSettings.log2e_start; log2eps <= CVGlobalSettings.log2e_end; log2eps = log2eps + CVGlobalSettings.log2e_incr) {
				
				// grid search for gamma
				for (int log2g = CVGlobalSettings.log2g_start; log2g <= CVGlobalSettings.log2g_end; log2g = log2g + CVGlobalSettings.log2g_incr) {
 					
					if (!silent) {
						System.out.println("C = 2^" + log2c + ", epsilon = 2^" + log2eps + ", gamma = 2^" + log2g);
					}
					
					String[] args2 = new String[15];
					args2[0] = "-v";
					args2[1] = ""+ CVGlobalSettings.folds;
					args2[2] = "-c";
					args2[3] = "" + Math.pow(2, log2c);
					args2[4] = "-s";
					args2[5] = "3";
					args2[6] = "-e";
					args2[7] = "" + Math.pow(2, log2eps);
					args2[8] = "-t";
					args2[9] = "2";
					args2[10] = "-g";
					args2[11] = "" + Math.pow(2, log2g);
					args2[12] = "-m";
					args2[13] = CVGlobalSettings.svm_cache + "";
					
					File f = new File(str_dir.trim() + "/" + str_file.trim());
					
					try {
						args2[14] = f.getCanonicalPath();
					} 
					catch (IOException e2) {
						e2.printStackTrace();
					}

					try {
						System.out.print("Command:   ./svm-train ");
						for (int t = 0; t < args2.length; t++) {
							System.out.print(args2[t] + " ");
						}
						System.out.println();

						/**
						 * evaluate over 10 randomly permutated cross-validation
						 * runs
						 */
						for (int runs = 0; runs < CVGlobalSettings.runs; runs++) {
							SeedFactory.setSeed(runs);
							svm_train.train(args2);
						}
						CVRegressionContainer cvrc = svm_train.getCvregcontainer();
						svm_train.clearCVResults();
						double[] all_q2 = cvrc.getAllQ2Values();
						double[] all_mse = cvrc.getAllMSEValues();
						double[] all_aae = cvrc.getAllAAEValues();

						/**
						 * compute means of mse, aae, q^2 and their stdv
						 */
						Mean mean = new Mean();
						double mse = mean.evaluate(all_mse);
						double q2 = mean.evaluate(all_q2);
						double aae = mean.evaluate(all_aae);
						StandardDeviation stdv = new StandardDeviation();
						double stdv_mse = stdv.evaluate(all_mse);
						double stdv_q2 = stdv.evaluate(all_q2);
						double stdv_aae = stdv.evaluate(all_aae); 

						
						/**
						 * append results
						 */
						if (!Double.isNaN(mse) && !Double.isNaN(q2)
								&& !Double.isNaN(stdv_q2)
								&& !Double.isNaN(stdv_q2)) {
							System.out.println("File:      " + cvrc.getMatrix());
							fw.append(cvrc.getMatrix() + "\t");
							
							System.out.println("C:         " + cvrc.getC());
							fw.append(cvrc.getC() + "\t");
							
							System.out.println("Gamma:     " + Math.pow(2, log2g));
							fw.append(Math.pow(2, log2g) + "\t");
							
							System.out.println("Epsilon:   " + cvrc.getEpsilon());
							fw.append(cvrc.getEpsilon() + "\t");
							
							System.out.println("MSE:       " + mse);
							fw.append(mse + "\t");
							
							System.out.println("MSE_STDV:  " + stdv_mse);
							fw.append(stdv_mse + "\t");
							
							System.out.println("AAE:       " + aae);
							fw.append(aae + "\t");
							
							System.out.println("AAE_STDV:  " + stdv_aae);
							fw.append(stdv_aae + "\t");
							
							System.out.println("Q2:        " + q2);
							fw.append(q2 + "\t");
							
							System.out.println("Q2_STDV:   " + stdv_q2 + "\n");
							fw.append(stdv_q2 + "\n");
							
							fw.flush();
						}

					} catch (IOException e1) {
						e1.printStackTrace();
						System.exit(1);
					}
				}
			}
		}
	}
}
