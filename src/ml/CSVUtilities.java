package ml;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.csv.CSVRecord;

public class CSVUtilities {
	
	/**
	 * Reads a csv whose first row is headers
	 * @param SAMPLE_CSV_FILE_PATH
	 * @param num_columns
	 * @param num_rows number of rows AFTER header
	 * @param column_names
	 * @return
	 * @throws IOException
	 */
	public static String[][] csv_reader_with_header(String SAMPLE_CSV_FILE_PATH,
			int num_columns, int num_rows, String[] column_names) throws IOException{
		
		Reader reader = Files.newBufferedReader(Paths.get(SAMPLE_CSV_FILE_PATH));
	    CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT
	    		.withFirstRecordAsHeader()
	            .withIgnoreHeaderCase()
	            .withTrim());
	    Iterable<CSVRecord> csvRecords = csvParser.getRecords();
	    String[][] output = new String [num_columns][num_rows];
	    int index = 0;
	    for (CSVRecord csvRecord : csvRecords) {
	    	for (int i=0; i < column_names.length; i++) {
	    		output[i][index] = csvRecord.get(column_names[i]);}
	    	index += 1;
	    }
	    csvParser.close();
		return transposer(output);
	}
	/**
	 * Reads a csv file without header
	 * @param SAMPLE_CSV_FILE_PATH
	 * @param num_columns
	 * @param num_rows
	 * @return
	 * @throws IOException
	 */
	public static String[][] csv_reader_without_header(String SAMPLE_CSV_FILE_PATH,
			int num_columns, int num_rows) throws IOException{
		
		Reader reader = Files.newBufferedReader(Paths.get(SAMPLE_CSV_FILE_PATH));
	    CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT
	    		.withFirstRecordAsHeader()
	            .withIgnoreHeaderCase()
	            .withTrim());
	    Iterable<CSVRecord> csvRecords = csvParser.getRecords();
	    String[][] output = new String [num_columns][num_rows];
	    int index = 0;
	    for (CSVRecord csvRecord : csvRecords) {
	    	for (int i=0; i < num_columns; i++) {
	    		output[i][index] = csvRecord.get(i);
	    	}
	    	index += 1;
	    }
	    csvParser.close();
		return transposer(output);
	}
	/**
	 * Transposes given string matrix
	 * @param in
	 * @return
	 */
	private static String[][] transposer(String[][] in){
		int numRows = in[0].length;
		int numCols = in.length;
		String[][] result = new String[numRows][numCols];
		for (int i = 0; i < numRows; i++) {
			for (int j = 0; j < numCols; j++) {
				result[i][j] = in[j][i];
			}
		}
		return result;
	}
	
	/**
	 Writes a csv file
	 @param {String} location - the folder directory in which the data will be saved
	 @param {String} name - name of the csv file that will be saved. Include '.csv' in the name
	 @param {String[]} column_names - Header names
	 @param {String[][]} data - the data matrix that will be converted to csv
	 @return {csv} - a csv file at the designated location
	 */
	public static void csv_writer_with_header(String location, String name, String[] column_names, String[][] data) throws IOException {
		
		FileWriter fileWriter = null;
		final String NEW_LINE_SEPARATOR = "\n";
		 CSVPrinter csvFilePrinter = null;
		 String fileName = location + name;
	
		//Create the CSVFormat object with "\n" as a record delimiter
		CSVFormat csvFileFormat = CSVFormat.DEFAULT.withRecordSeparator(NEW_LINE_SEPARATOR);
				
		//initialize FileWriter object
		fileWriter = new FileWriter(fileName);
		
		//initialize CSVPrinter object 
       csvFilePrinter = new CSVPrinter(fileWriter, csvFileFormat);
       
       //Create CSV file header
       List header = new ArrayList();
       for (int j=0; j < column_names.length; j++) {
       	header.add(column_names[j]);
       }
       csvFilePrinter.printRecord(header);
		
		//Write a new student object list to the CSV file
		for (int i = 0; i < data.length; i++) {
			List line = new ArrayList();
			for (int j=0; j < data[1].length; j++) {
				line.add(data[i][j]);
			}
           csvFilePrinter.printRecord(line);
		}
		csvFilePrinter.close();
	}
	/**
	 Writes a csv file
	 @param {String} location - the folder directory in which the data will be saved
	 @param {String} name - name of the csv file that will be saved. Include '.csv' in the name
	 @param {String[]} column_names - Header names
	 @param {String[][]} data - the data matrix that will be converted to csv
	 @return {csv} - a csv file at the designated location
	 */
	public static void csv_writer_with_header(String location, String name, String[] column_names, String[] data) throws IOException {
		
		FileWriter fileWriter = null;
		final String NEW_LINE_SEPARATOR = "\n";
		 CSVPrinter csvFilePrinter = null;
		 String fileName = location + name;
	
		//Create the CSVFormat object with "\n" as a record delimiter
		CSVFormat csvFileFormat = CSVFormat.DEFAULT.withRecordSeparator(NEW_LINE_SEPARATOR);
				
		//initialize FileWriter object
		fileWriter = new FileWriter(fileName);
		
		//initialize CSVPrinter object 
       csvFilePrinter = new CSVPrinter(fileWriter, csvFileFormat);
       
       //Create CSV file header
       List header = new ArrayList();
       for (int j=0; j < column_names.length; j++) {
       	header.add(column_names[j]);
       }
       csvFilePrinter.printRecord(header);
		
		//Write a new student object list to the CSV file
		for (int i = 0; i < data.length; i++) {
			List line = new ArrayList();
				line.add(data[i]);
           csvFilePrinter.printRecord(line);
		}
		csvFilePrinter.close();
	}
	/**
	 Writes a csv file
	 @param {String} location - the folder directory in which the data will be saved
	 @param {String} name - name of the csv file that will be saved. Include '.csv' in the name
	 @param {String[]} column_names - Header names
	 @param {double[][]} data - the data matrix that will be converted to csv
	 @return {csv} - a csv file at the designated location
	 */
	public static void csv_writer_with_header(String location, String name, String[] column_names, double[][] data) throws IOException {
		
		FileWriter fileWriter = null;
		final String NEW_LINE_SEPARATOR = "\n";
		 CSVPrinter csvFilePrinter = null;
		 String fileName = location + name;
	
		//Create the CSVFormat object with "\n" as a record delimiter
		CSVFormat csvFileFormat = CSVFormat.DEFAULT.withRecordSeparator(NEW_LINE_SEPARATOR);
				
		//initialize FileWriter object
		fileWriter = new FileWriter(fileName);
		
		//initialize CSVPrinter object 
       csvFilePrinter = new CSVPrinter(fileWriter, csvFileFormat);
       
       //Create CSV file header
       List header = new ArrayList();
       for (int j=0; j < column_names.length; j++) {
       	header.add(column_names[j]);
       }
       csvFilePrinter.printRecord(header);
		
		//Write a new student object list to the CSV file
		for (int i = 0; i < data.length; i++) {
			List line = new ArrayList();
			for (int j=0; j < data[0].length; j++) {
				line.add(data[i][j]);
			}
           csvFilePrinter.printRecord(line);
		}
		csvFilePrinter.close();
	}
}
