package ml;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
/**
 * Contains basic statistical methods like variance
 * @author K.Ataman
 *
 */
public class Statistics {
	public static double mean(double [] x) {
		double result = 0.0;
		for (int i = 0; i < x.length; i++) {
			result += x[i];
		}
		return result/x.length;
	}
	
	public static double variance(double[] x) {
		double mean = mean(x);
		double [] temp = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			temp[i] = (x[i] - mean) * (x[i] - mean);
		}
		return mean(temp);
	}
	
	public static double stdev(double[] x) {
		return Math.sqrt(variance(x));
	}
	
	public static double covariance(double[] x, double[] y){
		//error checking
		if (x.length==0) {
			throw new IllegalArgumentException("x is empty");
		}
		if (y.length==0) {
			throw new IllegalArgumentException("y is empty");
		}
		if (x.length!=y.length) {
			throw new IllegalArgumentException("Length of x and y must match");
		}
		
		double meanX = mean(x);
		double meanY = mean(y);
		double result = 0.0;
		for (int i = 0; i < x.length; i++) {
			result += (x[i] - meanX) * (y[i] - meanY);
		}
		return result / x.length;
	}
	
	public static double[][] covarianceMatrix(double[][]x){
		double[][] result = new double[x.length][x[0].length];
		for (int i = 0; i < x.length; i++) {
			for (int j = 0; j < x[0].length; j++) {
				result[i][j] = Statistics.covariance(x[i], x[j]);
			}
		}
		return result;
	}
	
	public static double pearsonCorrelationCoefficient(double[] x, double [] y){
		return Statistics.covariance(x, y) / (Statistics.stdev(x) * Statistics.stdev(y));
	}
	
	/**
	 * creates a pmf given an occurence array
	 * @param occurences all elements are nonnegative, is an array that has the occurences of 
	 * elements in it
	 * @param numDistinctElements number of different values occurencesarray can have
	 * @return pmf where pmf[i] is the chance of ith occurence happening
	 */
	public static double[] pmf(int[] occurences, int numDistinctElements) {
		double[] result = new double[numDistinctElements];
		//If an element occurs, increase the counter of that element in the results array
		for (int i = 0; i < occurences.length; i++) {
			result[occurences[i]]++;
		}
		
		//Normalize by dividing by the amount of occurences (instances)
		for (int i = 0; i < numDistinctElements; i++) {
			result[i] = (double) result[i] / occurences.length;
		}
		return result;
	}
	/**
	 * @param occurences format is [entry][variable]. So occurences[2][3] is the 3rd entry's 4th
	 * feature value
	 * @return a map where int[] is the elements' features' values, and the entry is the probability
	 * of that feature happening
	 */
	public static Map<int[], Double> jointPmf(int[][] occurences){
		
		//create and populate the instance array
		Map<int[], Double> pmf = new HashMap<int[], Double>();
		int[] key = new int[occurences[0].length];
		boolean keyFound = false;
		for (int i = 0; i < occurences.length; i++) {
			key = occurences[i];
			//might need .replace
			keyFound = pmf.containsKey(key);
			if (keyFound == true) {
				//the entry exists
				pmf.replace(key, pmf.get(key) + 1);
			}else {
				pmf.put(key, 1.0);
			}
		}
		
		//normalize
		for (Entry<int[], Double> entry: pmf.entrySet()) {
			entry.setValue((double)entry.getValue() / (double)occurences.length);
		}
		return pmf;
	}
}
