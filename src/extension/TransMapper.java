package extension;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.HashMap;

import resources.Resource;


public class TransMapper {
	
	public String getSuperclass(String key) {
		
		String infile = Resource.class.getResource("txt/map.ser").getFile();
		
		String res = "0.0.0.0.0";
		
		// load hm from file
		HashMap<String,String> hm = null;
		

		try {
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(infile));
			hm = (HashMap<String, String>) in.readObject();
			in.close();
		}
		catch(IOException ex) {
			ex.printStackTrace();
		}
		catch(ClassNotFoundException ex) {
			ex.printStackTrace();
		}

		// matching...
		for (int i=0; i<3; i++) {
			CharSequence key_tmp = key.toUpperCase().subSequence(0, key.length()-i);
			if (hm.containsKey(key_tmp)) {
				res = "" + hm.get(key_tmp);
				break;
			}
		}
		

		return res;
	}
	

	
	public static void main(String[] args) throws Exception {
		
		TransMapper mapper = new TransMapper();
		System.out.println(mapper.getSuperclass("CR"));
	}
}
