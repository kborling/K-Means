import java.io.File;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;
import java.util.ListIterator;

/*
 * @author Kevin Borling
 * CSCD 429 | Data Mining
 * Homework 2: Prediction of gene/protein localization
 * KDD Cup 2001 Challenge
 */
public class LocalizationMain 
{
	public static void main(String[] args) 
	{
		Scanner fin = null;
		LinkedList test = new LinkedList();
		LinkedList data = new LinkedList();
		String geneID, essential, classify, complex, phenotype, motif, chromosome, function, localization, line = "";

		try {

			// -------------- Read in Data Set -------------------//
			fin = new Scanner(new File("Genes_relation.data"));
			

			while(fin.hasNextLine()) 
			{
				// Remove trailing period
				line = fin.nextLine().replace(".","");
				// Parse File into tokens, splitting by commas or double quotations
	    		String [] tokens = line.split(",(?=([^\"]*\"[^\"]*\")*[^\"]*$)");
	    		geneID = tokens[0];
	    		essential = tokens[1];
	    		classify = tokens[2];
	    		complex = tokens[3];
	    		phenotype = tokens[4];
	    		motif = tokens[5];
	    		chromosome = tokens[6];
	    		function = tokens[7];
	    		localization = tokens[8];
	    		// Add values to data linked list
	    		data.add(geneID, essential, classify, complex, phenotype, motif, chromosome, function, localization);

			}// End while

			// --------------- Read in Test Set ------------------//
			fin = new Scanner(new File("Genes_relation.test"));

			while(fin.hasNextLine()) 
			{
				// Remove trailing period
				line = fin.nextLine().replace(".","");
				// Parse File into tokens, splitting by commas or double quotations
	    		String [] tokens = line.split(",(?=([^\"]*\"[^\"]*\")*[^\"]*$)");

	    		geneID = tokens[0];
	    		essential = tokens[1];
	    		classify = tokens[2];
	    		complex = tokens[3];
	    		phenotype = tokens[4];
	    		motif = tokens[5];
	    		chromosome = tokens[6];
	    		function = tokens[7];
	    		localization = tokens[8];
	    		// Add values to test linked list
	    		test.add(geneID, essential, classify, complex, phenotype, motif, chromosome, function, localization);
	    		
			}// End while

			// -------------- Sort elements in the list ---------------//
			test.sort();

			// ------------- Calculate Euclidean Distance -------------//
			test.getEuclid(data, test);

			// ------------- Remove duplicates -------------//
			test.removeDuplicates();

			// ------------- Print Results ---------------//
			System.out.print(test.toString());

		} catch(Exception e) {
		e.printStackTrace();
		}
				
	}// End main

}// End class