package ml;

import java.io.IOException;
//TODO add Bhattacharyya and mahanabolis distances
public class Distance {
	
	
	/**
	 * Checks if sizes of x and y are valid for distance calculation
	 * @param x {double[]}
	 * @param y {double[]}
	 * @throws IOException if lengths of x and y don't match or if x or y is empty
	 */
	private static void errorCheck(double [] x, double [] y){
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
	}
	/**
	 * Minkowski not included due to difficulties it presents
	 * @author K.Ataman
	 *
	 */
	public enum distanceMetric{
		EUCLIDIAN,
		CITY_BLOCK,
		//MINKOWSKI,
		//Mahanabolis
		SUP,
		PEARSON_CORRELATION
	}
	
	/**
	 * Finds the euclidian distance between x and y
	 * @param x : {double[]} the first element 
	 * @param y : {double[]} the second element 
	 * @return {double}
	 * @throws IOException if lenghts of x and y are erroneous
	 */
	public static double euclidian(double [] x, double [] y){
		return minkowski(2, x, y);
	}
	/**
	 * Finds the city block distance between x and y
	 * @param x : {double[]} the first element 
	 * @param y : {double[]} the second element 
	 * @return {double}
	 * @throws IOException if lenghts of x and y are erroneous
	 */
	public static double cityBlock(double [] x, double [] y){
		return minkowski(1, x, y);
	}
	/**
	 * Finds the minkowski distance between x and y
	 * @param power : {int} the power for which each individual distance will be raised to. For example, 
	 * if power = 2, then this becomes the euclidian distance
	 * @param x : {double[]} the first element 
	 * @param y : {double[]} the second element 
	 * @return {double}
	 * @throws IOException if lenghts of x and y are erroneous
	 */
	public static double minkowski(int power, double [] x, double [] y){
		//error checking
		errorCheck(x, y);
		
		double result = 0.0;
		for (int i=0; i < x.length; i++) {
			result += Math.pow(Math.abs(x[i] - y[i]), power);
		}
		return Math.pow(result, 1.0/power);
	}
	/**
	 * Finds the sup distance (max distance between features)
	 * @param x : {double[]} the first element 
	 * @param y : {double[]} the second element 
	 * @return {double}
	 * @throws IOException if lenghts of x and y are erroneous
	 */
	public static double sup(double [] x, double [] y){
		//error cehcking
		errorCheck(x, y);
		
		double max = 0.0;
		double temp = 0.0;
		for (int i=0; i < x.length; i++) {
			temp = Math.abs(x[i] - y[i]);
			if (temp > max) {
				max = temp;
			}
		}
		return max;
	}
	/**
	 * Based off on pearson correlation coefficient. Can detect dissimilarity, 
	 * but not the magnitude of differences. This is not a distance metric
	 * @param x {double []}
	 * @param y {double []}
	 * @return
	 * @throws IOException if lenghts of x and y are erroneous
	 */
	public static double pearsonCorrelation(double[] x, double[] y){
		//error checking
		errorCheck(x,y);

		return (1- Statistics.pearsonCorrelationCoefficient(x,y) ) / 2;
	}


/**
 * Calculates distance between x and y using the given distance metric
 * @param x
 * @param y
 * @param distanceMetric
 * @return
 * @throws IOException
 */
public static double calculateDistance(double[] x, double[] y, Distance.distanceMetric distanceMetric
		){
	switch(distanceMetric) {
		case EUCLIDIAN:
			return Distance.euclidian(x, y);
		case CITY_BLOCK:
			return Distance.cityBlock(x, y);
		case SUP:
			return Distance.sup(x, y);
		case PEARSON_CORRELATION:
			return Distance.pearsonCorrelation(x, y);
	}
	return 0;
	}
}
