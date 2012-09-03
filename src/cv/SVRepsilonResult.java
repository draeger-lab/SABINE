package cv;

public class SVRepsilonResult {

	private String str_matrix;
	private double q2;
	private double q2_stdv;
	private double rmse;
	private double rmse_stdv;
	private double c;
	private double epsilon;

	/**
	 * check wether one of the results is NaN or Infinity
	 * 
	 * @return
	 */
	public boolean checkIntegrity() {
		if (checkInfinitity(q2) || checkInfinitity(q2_stdv) || checkInfinitity(rmse) || checkInfinitity(rmse_stdv)
				|| checkNaN(q2) || checkNaN(q2_stdv) || checkNaN(rmse) || checkNaN(rmse_stdv)) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * checks wether a double is inifinity
	 * 
	 * @param d
	 * @return
	 */
	private boolean checkInfinitity(double d) {
		return Double.isInfinite(d);
	}

	/**
	 * checks wether a double is NaN
	 * 
	 * @param d
	 * @return
	 */
	private boolean checkNaN(double d) {
		return Double.isNaN(d);
	}

	public double getC() {
		return c;
	}

	public void setC(double c) {
		this.c = c;
	}

	public double getEpsilon() {
		return epsilon;
	}

	public void setEpsilon(double epsilon) {
		this.epsilon = epsilon;
	}

	public double getQ2() {
		return q2;
	}

	public void setQ2(double q2) {
		this.q2 = q2;
	}

	public double getQ2_stdv() {
		return q2_stdv;
	}

	public void setQ2_stdv(double q2_stdv) {
		this.q2_stdv = q2_stdv;
	}

	public double getRmse() {
		return rmse;
	}

	public void setRmse(double rmse) {
		this.rmse = rmse;
	}

	public double getRmse_stdv() {
		return rmse_stdv;
	}

	public void setRmse_stdv(double rmse_stdv) {
		this.rmse_stdv = rmse_stdv;
	}

	public String getStr_matrix() {
		return str_matrix;
	}

	public void setStr_matrix(String str_matrix) {
		this.str_matrix = str_matrix;
	}

}
