package cv;

import java.util.Random;

public class SeedFactory {

	/**
	 * Singleton class for distribution of equal seeds
	 */
	private static SeedFactory sf;
	private static Random generator = new Random();

	private SeedFactory() {
	}

	public static void setSeed(int init) {
		generator.setSeed(init);
	}

	public static double getNewSeededRandomNumber() {
		return generator.nextDouble() / Double.MAX_VALUE;
	}

	public static SeedFactory getInstance() {
		if (sf == null) {
			sf = new SeedFactory();
			return sf;
		} else {
			return sf;
		}
	}
}
