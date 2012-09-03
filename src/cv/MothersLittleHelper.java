package cv;

public class MothersLittleHelper {

	public static Double[] doubleArray2DoubleArray(double[] d) {
		Double[] result = new Double[d.length];
		for (int i = 0; i < d.length; i++) {
			result[i] = new Double(d[i]);
		}
		return result;
	}

	public static double[] DoubleArray2doubleArray(Double[] d) {
		double[] result = new double[d.length];
		for (int i = 0; i < d.length; i++) {
			result[i] = d[i].doubleValue();
		}
		return result;
	}

}
