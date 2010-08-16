package util.parsers.tools.IPRscan;
import java.util.ArrayList;


public class MainIPRscan {
	
	
public static void main(String[] args) {
		
		//// Input
		String infile = args[0];
		String outfile = args[1];

		InputFileParser inputfileparser = new InputFileParser();
		inputfileparser.parse(infile);
		
		ArrayList<String> symbols = inputfileparser.getSymbols();
		ArrayList<ArrayList<String>> uniprot_ids = inputfileparser.getUniprot_ids();
		ArrayList<String> species = inputfileparser.getSpecies();
		ArrayList<ArrayList<ArrayList<String>>> matrices = inputfileparser.getMatrices();
		ArrayList<ArrayList<String>> sequences = inputfileparser.getSequences();
		ArrayList<ArrayList<String>> bindings = inputfileparser.getBindings();
		ArrayList<String> transfac = inputfileparser.getTransfac();
		ArrayList<Boolean> flags_arr = inputfileparser.getFlags_arr();
		ArrayList<ArrayList<String>> MN = inputfileparser.getMN();
		
		
		////
		Converter converter = new Converter();
		// convert matrices from absolute to relative (depreciated)
		//matrices = converter.a2r(matrices);
		
		//// Filter
		// filter remaining flagged entries
		Filter filter = new Filter();
		filter.checkFlag(symbols, uniprot_ids, species, matrices, sequences, bindings, transfac, flags_arr, MN);
		
		//// Output
		// convert ArrayList to Array
		boolean[] flags = converter.l2r(flags_arr);		
		
		// write
		OutputFileWriter outputfilewriter = new OutputFileWriter();
		outputfilewriter.writeOutfile(symbols, uniprot_ids, species, matrices, sequences, bindings, transfac, flags, MN, outfile);

	}
}


// TODO:
// IPRscan for entries marked by "IP" -> (flag == true)
