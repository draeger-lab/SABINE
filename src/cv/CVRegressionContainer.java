package cv;

import java.util.ArrayList;

public class CVRegressionContainer {

	private String matrix;
	private int folds;
	private ArrayList<Double[]> q2;
	private ArrayList<Double[]> mse;
	private ArrayList<Double[]> aae;
	private boolean integrity = true;
	private double C;
	private double epsilon;

	public int getFolds() {
		return folds;
	}

	public void setFolds(int folds) {
		this.folds = folds;
	}

	public void addToQ2(Double[] d) {
		this.q2.add(d);
	}

	public void addToMSE(Double[] d) {
		this.mse.add(d);
	}

	public void addToAAE(Double[] d){
		this.aae.add(d);
	}
	
	public CVRegressionContainer() {
		init();
	}

	public ArrayList<Double[]> getQ2() {
		return q2;
	}

	public String getMatrix() {
		return matrix;
	}

	public void setMatrix(String matrix) {
		this.matrix = matrix;
	}

	public void init() {
		this.q2 = new ArrayList<Double[]>();
		this.mse = new ArrayList<Double[]>();
		this.aae = new ArrayList<Double[]>();
		this.integrity = true;
	}

	/**
	 * 
	 * gets all squared corellation coefficients (Q^2)
	 * 
	 * @return
	 */
	public double[] getAllQ2Values() {
		ArrayList<Double> result = new ArrayList<Double>();
		for (int i = 0; i < this.q2.size(); i++) {
			Double[] al = this.q2.get(i);
			for (int j = 0; j < al.length; j++) {
				result.add(al[j]);
			}
		}
		double[] r = new double[result.size()];
		for (int i = 0; i < result.size(); i++) {
			r[i] = result.get(i).doubleValue();
		}
		return r;
	}

	/**
	 * gets all squared errors (MSE)
	 * 
	 * @return
	 */
	public double[] getAllMSEValues() {
		ArrayList<Double> result = new ArrayList<Double>();
		for (int i = 0; i < this.mse.size(); i++) {
			Double[] al = this.mse.get(i);
			for (int j = 0; j < al.length; j++) {
				result.add(al[j]);
			}
		}
		double[] r = new double[result.size()];
		for (int i = 0; i < result.size(); i++) {
			r[i] = result.get(i).doubleValue();
		}
		return r;
	}
	
	/**
	 * gets all absolute errors
	 * 
	 * @return
	 */
	public double[] getAllAAEValues() {
		ArrayList<Double> result = new ArrayList<Double>();
		for (int i = 0; i < this.aae.size(); i++) {
			Double[] al = this.aae.get(i);
			for (int j = 0; j < al.length; j++) {
				result.add(al[j]);
			}
		}
		double[] r = new double[result.size()];
		for (int i = 0; i < result.size(); i++) {
			r[i] = result.get(i).doubleValue();
		}
		return r;
	}	

	public double getC() {
		return C;
	}

	public void setC(double c) {
		C = c;
	}

	public double getEpsilon() {
		return epsilon;
	}

	public void setEpsilon(double epsilon) {
		this.epsilon = epsilon;
	}

	public boolean isIntegrity() {
		return integrity;
	}
}
