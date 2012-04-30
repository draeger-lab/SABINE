package libsvm;

public class LabeledTF implements Comparable {
	
	String name = null;
	
	double label = 0.0;
	
	
	public LabeledTF(String n, double l) {

		name = n;
		label = l;
	
	}
	
	public String getName() {
		return this.name;
	}
	
	public double getLabel() {
		return this.label;
	}


	public int compareTo(Object o) {
		
		
		LabeledTF tf = (LabeledTF) o;
		
		
		if(label < tf.label) {
			
			return 1;
		
		}
		
		else {
			
			if(label > tf.label) {
				
				return -1;
				
			}
			
			else {
				
				return 0;
			}

		}
		
	}
	
}
