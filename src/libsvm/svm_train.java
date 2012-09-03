package libsvm;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.Vector;

import cv.CVGlobalSettings;
import cv.MothersLittleHelper;
import cv.CVRegressionContainer;

import libsvmcore.svm;
import libsvmcore.svm_model;
import libsvmcore.svm_node;
import libsvmcore.svm_parameter;
import libsvmcore.svm_problem;

public class svm_train {
	
	private boolean silent = CVGlobalSettings.silent;
	
	private svm_parameter param; // set by parse_command_line
	private svm_problem prob; // set by read_problem
	private svm_model model;
	/**
	 * 
	 */
	private String input_file_name; // set by parse_command_line
	private String model_file_name; // set by parse_command_line
	private String error_msg;
	private int cross_validation;
	private int nr_fold;
	private static CVRegressionContainer cvregcontainer = new CVRegressionContainer();
	private static ArrayList<Double> accuracy = new ArrayList<Double>();

	private static void exit_with_help() {
		System.out
				.print("Usage: svm_train [options] training_set_file [model_file]\n"
						+ "options:\n"
						+ "-s svm_type : set type of SVM (default 0)\n"
						+ "	0 -- C-SVC\n"
						+ "	1 -- nu-SVC\n"
						+ "	2 -- one-class SVM\n"
						+ "	3 -- epsilon-SVR\n"
						+ "	4 -- nu-SVR\n"
						+ "-t kernel_type : set type of kernel function (default 2)\n"
						+ "	0 -- linear: u'*v\n"
						+ "	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
						+ "	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
						+ "	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
						+ "	4 -- precomputed kernel (kernel values in training_set_file)\n"
						+ "-d degree : set degree in kernel function (default 3)\n"
						+ "-g gamma : set gamma in kernel function (default 1/k)\n"
						+ "-r coef0 : set coef0 in kernel function (default 0)\n"
						+ "-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
						+ "-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
						+ "-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
						+ "-m cachesize : set cache memory size in MB (default 100)\n"
						+ "-e epsilon : set tolerance of termination criterion (default 0.001)\n"
						+ "-h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
						+ "-b probability_estimates: whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
						+ "-wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
						+ "-v n: n-fold cross validation mode\n");
		System.exit(1);
	}

	private void do_cross_validation() {
		int total_correct = 0;
		double total_error = 0;
		double absolute_error = 0;
		double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
		double[] target = new double[prob.l];

		ArrayList<Double[]> folds = svm.svm_cross_validation(prob, param, nr_fold, target);

		if (param.svm_type == svm_parameter.EPSILON_SVR
				|| param.svm_type == svm_parameter.NU_SVR) {
			for (int i = 0; i < prob.l; i++) {
				// label
				double y = prob.y[i];
				// predictions with cross-validation
				double v = target[i];
				// absolute error of prediction
				total_error = total_error + Math.abs(v - y) * Math.abs(v - y);
				absolute_error = absolute_error + Math.abs(v - y);
				sumv += v;
				sumy += y;
				sumvv += v * v;
				sumyy += y * y;
				sumvy += v * y;
			}

			// Q^2, RMSE from k-folds
			// cvregcontainer.addToTotalError(total_error / prob.l);
			double[] foldwise_r2 = new double[nr_fold];
			double[] foldwise_mse = new double[nr_fold];
			double[] foldwise_aae = new double[nr_fold];

			for (int k = 0; k < folds.size(); k++) {

				Double[] fold_res = folds.get(k);
				// get label and predictions from this fold
				Vector<Double> pred = new Vector<Double>();
				Vector<Double> targ = new Vector<Double>();
				
				
				if(!silent)	System.out.println("Y(=label)\tX(=prediction)\tfold=" + (k + 1));
				
				
				for (int p = 0; p < fold_res.length; p++) {
					if (fold_res[p] != null) {
						
						if(!silent)	System.out.println(CVGlobalSettings.df.format(prob.y[p]) + "\t" + CVGlobalSettings.df.format(fold_res[p]));
						
						double tar = prob.y[p];
						double pre = fold_res[p];
						pred.add(pre);
						targ.add(tar);
					}
				}
				foldwise_r2[k] = getQ24Fold(targ, pred);
				
				
				if(!silent)	System.out.println("\nSUMMARY(fold=" + (k + 1) + ")");
				if(!silent)	System.out.println("R^2\t" + CVGlobalSettings.df.format(foldwise_r2[k]));
				
				
				foldwise_mse[k] = getMSE4Fold(targ, pred);
				
				if(!silent)	System.out.println("MSE\t" + CVGlobalSettings.df.format(foldwise_mse[k]));
				
				foldwise_aae[k] = getAAEFold(targ, pred);
				
				if(!silent)	System.out.println("AAE\t" + CVGlobalSettings.df.format(foldwise_aae[k]) + "\n");

			}

			cvregcontainer.addToMSE(MothersLittleHelper.doubleArray2DoubleArray(foldwise_mse));
			cvregcontainer.addToQ2(MothersLittleHelper.doubleArray2DoubleArray(foldwise_r2));
			cvregcontainer.addToAAE(MothersLittleHelper.doubleArray2DoubleArray(foldwise_aae));
			cvregcontainer.setC(param.C);
			cvregcontainer.setEpsilon(param.eps);
			cvregcontainer.setFolds(nr_fold);
			cvregcontainer.setMatrix(input_file_name);

		} else {
			for (int i = 0; i < prob.l; i++)
				if (target[i] == prob.y[i]) {
					++total_correct;
				}
			/**
			 * ACCURACY
			 */
			double acc = 100.0 * total_correct / prob.l;
			System.out.print("Cross Validation Accuracy = " + acc + "%\n");
			/**
			 * TODO SENSITIVITY, SPECIFITY
			 */
			svm_train.accuracy.add(acc);
		}
	}

	/**
	 * clears the accuracy
	 */
	public static void clearAccuracies() {
		accuracy = new ArrayList<Double>();
	}

	/**
	 * 
	 * computes q^2
	 * 
	 * @param targ
	 * @param pred
	 * @return
	 */
	private static double getQ24Fold(Vector<Double> targ, Vector<Double> pred) {
		// now compute the r2 from the kth fold
		double total_error = 0;
		double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
		for (int i = 0; i < targ.size(); i++) {
			// label
			double y = targ.get(i);
			// predictions with cross-validation
			double v = pred.get(i);
			// absolute error of prediction
			total_error += (v - y) * (v - y);
			sumv += v;
			sumy += y;
			sumvv += v * v;
			sumyy += y * y;
			sumvy += v * y;
		}
		double r2 = (targ.size() * sumvy - sumv * sumy)
				* (targ.size() * sumvy - sumv * sumy)
				/ ((targ.size() * sumvv - sumv * sumv) * (targ.size() * sumyy - sumy
						* sumy));
		return r2;
	}

	/**
	 * 
	 * computes RMSE for each fold
	 * 
	 * @param targ
	 * @param pred
	 * @return
	 */
	private static double getMSE4Fold(Vector<Double> targ, Vector<Double> pred) {
		double result = 0;
		for (int i = 0; i < targ.size(); i++) {
			result = result + Math.abs(targ.get(i) - pred.get(i))
					* Math.abs(targ.get(i) - pred.get(i));
		}
		return result / targ.size();
	}

	/**
	 * 
	 * computes AAE for each fold
	 * 
	 * @param targ
	 * @param pred
	 * @return
	 */
	private static double getAAEFold(Vector<Double> targ, Vector<Double> pred) {
		double result = 0;
		for (int i = 0; i < targ.size(); i++) {
			result = result + Math.abs(targ.get(i) - pred.get(i));
		}
		return result / targ.size();
	}

	private void run(String argv[]) throws IOException {
		parse_command_line(argv);
		read_problem();
		error_msg = svm.svm_check_parameter(prob, param);

		if (error_msg != null) {
			System.err.print("Error: " + error_msg + "\n");
			System.exit(1);
		}

		if (cross_validation != 0) {
			do_cross_validation();
		} else {
			model = svm.svm_train(prob, param);
			svm.svm_save_model(model_file_name, model);
		}
	}

	public static void main(String[] argv) throws IOException
	{
		svm_train t = new svm_train();
		
		t.run(argv);
	}
	
	public static CVRegressionContainer train(String argv[]) throws IOException {
		svm_train t = new svm_train();
		t.run(argv);
		return cvregcontainer;
	}

	private static double atof(String s) {
		return Double.valueOf(s).doubleValue();
	}

	private static int atoi(String s) {
		return Integer.parseInt(s);
	}

	public static ArrayList<Double> getAccuracy() {
		return accuracy;
	}

	private void parse_command_line(String argv[]) {
		int i;

		param = new svm_parameter();
		// default values
		param.svm_type = svm_parameter.C_SVC;
		param.kernel_type = svm_parameter.RBF;
		param.degree = 3;
		param.gamma = 0; // 1/k
		param.coef0 = 0;
		param.nu = 0.5;
		param.cache_size = 100;
		param.C = 1;
		param.eps = 1e-3;
		param.p = 0.1;
		param.shrinking = 1;
		param.probability = 0;
		param.nr_weight = 0;
		param.weight_label = new int[0];
		param.weight = new double[0];
		cross_validation = 0;

		// parse options
		for (i = 0; i < argv.length; i++) {
			if (argv[i].charAt(0) != '-')
				break;
			if (++i >= argv.length)
				exit_with_help();
			switch (argv[i - 1].charAt(1)) {
			case 's':
				param.svm_type = atoi(argv[i]);
				break;
			case 't':
				param.kernel_type = atoi(argv[i]);
				break;
			case 'd':
				param.degree = atoi(argv[i]);
				break;
			case 'g':
				param.gamma = atof(argv[i]);
				break;
			case 'r':
				param.coef0 = atof(argv[i]);
				break;
			case 'n':
				param.nu = atof(argv[i]);
				break;
			case 'm':
				param.cache_size = atof(argv[i]);
				break;
			case 'c':
				param.C = atof(argv[i]);
				break;
			case 'e':
				param.eps = atof(argv[i]);
				break;
			case 'p':
				param.p = atof(argv[i]);
				break;
			case 'h':
				param.shrinking = atoi(argv[i]);
				break;
			case 'b':
				param.probability = atoi(argv[i]);
				break;
			case 'v':
				cross_validation = 1;
				nr_fold = atoi(argv[i]);
				if (nr_fold < 2) {
					System.err.print("n-fold cross validation: n must >= 2\n");
					exit_with_help();
				}
				break;
			case 'w':
				++param.nr_weight;
				{
					int[] old = param.weight_label;
					param.weight_label = new int[param.nr_weight];
					System.arraycopy(old, 0, param.weight_label, 0,
							param.nr_weight - 1);
				}

				{
					double[] old = param.weight;
					param.weight = new double[param.nr_weight];
					System.arraycopy(old, 0, param.weight, 0,
							param.nr_weight - 1);
				}

				param.weight_label[param.nr_weight - 1] = atoi(argv[i - 1]
						.substring(2));
				param.weight[param.nr_weight - 1] = atof(argv[i]);
				break;
			default:
				System.err.print("unknown option\n");
				exit_with_help();
			}
		}

		// determine filenames

		if (i >= argv.length)
			exit_with_help();

		input_file_name = argv[i];

		if (i < argv.length - 1)
			model_file_name = argv[i + 1];
		else {
			int p = argv[i].lastIndexOf('/');
			++p; // whew...
			model_file_name = argv[i].substring(p) + ".model";
		}
	}

	// read in a problem (in svmlight format)

	private void read_problem() throws IOException {

		BufferedReader fp = new BufferedReader(new FileReader(input_file_name));
		Vector<String> vy = new Vector<String>();
		Vector<svm_node[]> vx = new Vector<svm_node[]>();
		int max_index = 0;

		while (true) {
			String line = fp.readLine();
			if (line == null)
				break;

			StringTokenizer st = new StringTokenizer(line, " \t\n\r\f:");

			vy.addElement(st.nextToken());
			int m = st.countTokens() / 2;
			svm_node[] x = new svm_node[m];
			for (int j = 0; j < m; j++) {
				x[j] = new svm_node();
				x[j].index = atoi(st.nextToken());
				x[j].value = atof(st.nextToken());
			}
			if (m > 0)
				max_index = Math.max(max_index, x[m - 1].index);
			vx.addElement(x);
		}

		prob = new svm_problem();
		prob.l = vy.size();
		prob.x = new svm_node[prob.l][];
		for (int i = 0; i < prob.l; i++)
			prob.x[i] = vx.elementAt(i);
		prob.y = new double[prob.l];
		for (int i = 0; i < prob.l; i++)
			prob.y[i] = atof(vy.elementAt(i));

		if (param.gamma == 0)
			param.gamma = 1.0 / max_index;

		if (param.kernel_type == svm_parameter.PRECOMPUTED)
			for (int i = 0; i < prob.l; i++) {
				if (prob.x[i][0].index != 0) {
					System.err
							.print("Wrong kernel matrix: first column must be 0:sample_serial_number\n");
					System.exit(1);
				}
				if ((int) prob.x[i][0].value <= 0
						|| (int) prob.x[i][0].value > max_index) {
					System.err
							.print("Wrong input format: sample_serial_number out of range\n");
					System.exit(1);
				}
			}

		fp.close();
	}

	public static CVRegressionContainer getCvregcontainer() {
		return cvregcontainer;
	}

	public static void clearCVResults() {
		cvregcontainer = new CVRegressionContainer();
	}
}
