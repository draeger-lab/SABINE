package cv;

/**
 * 
 * Test class for seeded random number generator
 * 
 * @author hinselma
 *
 */
public class SeedFactoryTest {
	
	/**
	 * @param args
	 */
	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		SeedFactory sf;
		sf = SeedFactory.getInstance();

		sf.setSeed(1);
		for (int i = 0; i < 10; i++) {
			System.out.println(sf.getNewSeededRandomNumber());
		}
		System.out.println("");
		sf.setSeed(2);
		for (int i = 0; i < 10; i++) {
			System.out.println(sf.getNewSeededRandomNumber());
		}
		System.out.println("");
		sf.setSeed(1);
		for (int i = 0; i < 10; i++) {
			System.out.println(sf.getNewSeededRandomNumber());
		}
	}
}
