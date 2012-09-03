package cv;

public class GridSizeTester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CVGlobalSettings cvgrid = new CVGlobalSettings();
		System.out.println("\nParameter grids:");
		System.out.println("");		
		cvgrid.printCParameterGrid();
		cvgrid.printGammaParameterGrid();
		cvgrid.printEpsilonParameterGrid();
		cvgrid.printConfig();
		System.out.println("");
	//	System.out.println("");
	//	System.out.println("Classification (Precomputed):\t"+ cvgrid.getKernelGridSizeClassifaction());
	//	System.out.println("Regression (Precomputed):\t"+ cvgrid.getKernelGridSizeRegression());
	//	System.out.println("Classification (RBF):\t\t"+ cvgrid.getRBFGridSizeClassifaction());
		System.out.println("Total grid size:\t\t:\t" + cvgrid.getRBFGridSizeRegression() + "\n");
	}
}
