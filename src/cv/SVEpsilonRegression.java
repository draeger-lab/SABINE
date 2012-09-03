package cv;

import libsvm.svm_train;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;


public class SVEpsilonRegression {

	/**
	 * 
	 * epsilon regression
	 * 
	 * @param str_files
	 * @param str_dir
	 */
	public static void performEpsilonRegression(String str_file, String str_dir,
			FileWriter fw) {

		// grid search for C parameter
		for (int log2c = CVGlobalSettings.log2c_start; log2c <= CVGlobalSettings.log2c_end; log2c = log2c + CVGlobalSettings.log2c_incr) {
 
			// grid search for epsilon
			for (int log2eps = CVGlobalSettings.log2e_start; log2eps <= CVGlobalSettings.log2e_end; log2eps = log2eps + CVGlobalSettings.log2e_incr) {

 				String[] args2 = new String[13];
				args2[0] = "-v";
				args2[1] = ""+CVGlobalSettings.folds;
				args2[2] = "-c";
				args2[3] = "" + Math.pow(2, log2c);
				args2[4] = "-s";
				args2[5] = "3";
				args2[6] = "-e";
				args2[7] = "" + Math.pow(2, log2eps);
				args2[8] = "-t";
				args2[9] = "4";
				args2[10] = "-m";
				args2[11] = CVGlobalSettings.svm_cache + "";
				File f = new File(str_dir.trim() + "/" + str_file.trim());
				try {
					args2[12] = f.getCanonicalPath();
				} catch (IOException e2) {
					e2.printStackTrace();
				}

				try {
					for (int t = 0; t < args2.length; t++) {
						System.out.print(args2[t] + " ");
					}
					System.out.println("");

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
							&& !Double.isNaN(stdv_q2) && !Double.isNaN(stdv_q2)) {
						System.out.println("CV-SUMMARY, folds="+(CVGlobalSettings.folds)+", runs="+CVGlobalSettings.runs);
						System.out.println("File:\t\t" + cvrc.getMatrix());
						fw.append(cvrc.getMatrix() + "\t");
						System.out.println("log2C:\t\t" + log2c);
						fw.append(log2c + "\t");
						System.out.println("log2e:\t\t" + log2eps);
						fw.append(log2eps + "\t");
						System.out.println("MSE:\t\t" + CVGlobalSettings.df.format(mse));
						fw.append(CVGlobalSettings.df.format(mse) + "\t");
						System.out.println("MSE_STDV:\t" + CVGlobalSettings.df.format(stdv_mse));
						fw.append(CVGlobalSettings.df.format(stdv_mse) + "\t");
						System.out.println("AAE:\t\t" + CVGlobalSettings.df.format(aae));
						fw.append(CVGlobalSettings.df.format(aae) + "\t");
						System.out.println("AAE_STDV:\t" + CVGlobalSettings.df.format(stdv_aae));
						fw.append(CVGlobalSettings.df.format(stdv_aae) + "\t");
						System.out.println("Q2:\t\t" + CVGlobalSettings.df.format(q2));
						fw.append(CVGlobalSettings.df.format(q2) + "\t");
						System.out.println("Q2_STDV:\t" + CVGlobalSettings.df.format(stdv_q2));
						fw.append(CVGlobalSettings.df.format(stdv_q2) + "\n");
						System.out.println("");
					}

				} catch (IOException e1) {
					e1.printStackTrace();
					System.exit(1);
				}
			}
		}
	}
}
