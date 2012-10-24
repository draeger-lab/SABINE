package cv;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;

public class SVClassification {

	/**
	 * 
	 * epsilon regression
	 * 
	 * @param str_files
	 * @param str_dir
	 */
	public static void performClassification(String str_file, String str_dir,
			FileWriter fw) {

		// grid search for C parameter
		for (int log2c = CVGlobalSettings.log2c_start; log2c <= CVGlobalSettings.log2c_end; log2c = log2c + CVGlobalSettings.log2c_incr) {

			String[] args2 = new String[11];
			args2[0] = "-v";
			args2[1] = ""+CVGlobalSettings.folds;
			args2[2] = "-c";
			args2[3] = "" + Math.pow(2, log2c);
			args2[4] = "-s";
			args2[5] = "0";
			args2[6] = "-t";
			args2[7] = "4";
			args2[8] = "-m";
			args2[9] = CVGlobalSettings.svm_cache + "";
			File f = new File(str_dir.trim() + "/" + str_file.trim());
			try {
				args2[10] = f.getCanonicalPath();
			} catch (IOException e2) {
				e2.printStackTrace();
			}

			try {
				for (int t = 0; t < args2.length; t++) {
					System.out.print(args2[t] + " ");
				}
				System.out.println("");

				/**
				 * evaluate over 10 randomly permutated cross-validation runs
				 */
				for (int runs = 0; runs < CVGlobalSettings.runs; runs++) {
					SeedFactory.setSeed(runs);
					libsvm.svm_train.train(args2);
				}
				ArrayList<Double> result = libsvm.svm_train.getAccuracy();

				double[] accuracies = new double[result.size()];
				for (int r = 0; r < result.size(); r++) {
					accuracies[r] = result.get(r);
				}
				libsvm.svm_train.clearAccuracies();
				System.out.println("");

				/**
				 * compute accuracy and its stdv
				 */
				Mean mean = new Mean();
				double acc_mean = mean.evaluate(accuracies);
				StandardDeviation sdtv = new StandardDeviation();
				double acc_stdv = sdtv.evaluate(accuracies);

				/**
				 * append results
				 */
				if (!Double.isNaN(acc_mean) && !Double.isNaN(acc_stdv)) {
					System.out.println("FILE:\t\t" + f.getCanonicalPath());
					fw.append(f.getCanonicalPath() + "\t");
					System.out.println("log2c:\t\t" + log2c);
					fw.append(log2c + "\t");
					System.out.println("ACC MEAN:\t" + CVGlobalSettings.df.format(acc_mean));
					fw.append(CVGlobalSettings.df.format(acc_mean) + "\t");
					System.out.println("ACC STDV:\t" + CVGlobalSettings.df.format(acc_stdv));
					fw.append(CVGlobalSettings.df.format(acc_stdv) + "\n");
					System.out.println("");
				}

			} catch (IOException e1) {
				e1.printStackTrace();
				System.exit(1);
			}
		}
	}
}