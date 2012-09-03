package cv;

import java.text.DecimalFormat;

/**
 * 
 * Global access for all important parameters for grid search
 * 
 * @author hinselma
 * 
 */
public class CVGlobalSettings {
	
	
	public static boolean silent = true;
	public static boolean verbose = false;

	public static DecimalFormat df = new DecimalFormat("#.##");
	
	/**
	 * threshold of different labels. if a problem has more labels than this
	 * threshold, the problem will be treated as regression problem
	 */
	public static int regression_threshold = 5;

	public static int svm_cache = 500;

	/**
	 * set the grid for log2c (all kernels, classification and regression)
	 */
	public static int log2c_start = -4;
	public static int log2c_end = 10;
	public static int log2c_incr = 1;

	public void printCParameterGrid() {
		System.out.println("log2(C) \t Start\t\t:\t" + log2c_start);
		System.out.println("log2(C) \t End\t\t:\t"   + log2c_end);
		System.out.println("log2(C) \t Increment\t:\t"  + log2c_incr);
	}

	/**
	 * set range for log2gamma (RBF kernel implemented in LIBSVM, classification
	 * and regression)
	 */
	public static int log2g_start = -14;
	public static int log2g_end = -1;
	public static int log2g_incr = 2;

	public void printGammaParameterGrid() {
		System.out.println("\nlog2(gamma) \t Start\t\t:\t" + log2g_start);
		System.out.println("log2(gamma) \t End\t\t:\t" + log2g_end);
		System.out.println("log2(gamma) \t Increment\t:\t" + log2c_incr);
	}

	/**
	 * set range for log2epsilon for epsilon regression (epsilon regression)
	 */
	public static int log2e_start = -8;
	public static int log2e_end = -1;
	public static int log2e_incr = 1;

	public void printEpsilonParameterGrid() {
		System.out.println("\nlog2(epsilon) \t Start\t\t:\t" + log2e_start);
		System.out.println("log2(epsilon) \t End\t\t:\t" + log2e_end);
		System.out.println("log2(epsilon) \t Increment\t:\t" + log2e_incr);
	}

	/**
	 * runs and number of folds in a k-fold cross-validation run
	 */
	public static int runs = 1;
	public static int folds = 5;

	public void printConfig() {
		System.out.println("\n#runs\t\t\t\t:\t" + runs);
		System.out.println("#folds per run\t\t\t:\t" + folds);
	}

	public int getRBFGridSizeClassifaction() {
		int cgrid = Math.abs(log2c_end) + Math.abs(log2c_start);
		int ggrid = Math.abs(log2g_end) + Math.abs(log2g_start);
		return cgrid * ggrid / (log2c_incr * log2g_incr);
	}

	public int getRBFGridSizeRegression() {
		int cgrid = (log2c_end - log2c_start) / log2c_incr + 1;
		int egrid = (log2e_end - log2e_start) / log2e_incr + 1;
		int ggrid = (log2g_end - log2g_start) / log2g_incr + 1;
		return cgrid * egrid * ggrid;
	}

	public int getKernelGridSizeRegression() {
		int cgrid = Math.abs(log2c_end) + Math.abs(log2c_start);
		int egrid = Math.abs(log2e_end) + Math.abs(log2e_start);
		return cgrid * egrid / (log2c_incr * log2e_incr);
	}

	public int getKernelGridSizeClassifaction() {
		int cgrid = Math.abs(log2c_end) + Math.abs(log2c_start);
		return cgrid / (log2c_incr);
	}
}
