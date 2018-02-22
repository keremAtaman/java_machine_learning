package ml;

import java.io.IOException;

public class Distance {
	/**
	 * Checks if sizes of x and y are valid for distance calculation
	 * @param x {double[]}
	 * @param y {double[]}
	 * @throws IOException if lengths of x and y don't match or if x or y is empty
	 */
	private static void errorCheck(double [] x, double [] y) throws IOException{
		//error checking
		if (x.length==0) {
			throw new IOException("x is empty");
		}
		if (y.length==0) {
			throw new IOException("y is empty");
		}
		if (x.length!=y.length) {
			throw new IOException("Length of x and y must match");
		}
	}
	/**
	 * Finds the euclidian distance between x and y
	 * @param x : {double[]} the first element 
	 * @param y : {double[]} the second element 
	 * @return {double}
	 * @throws IOException if lenghts of x and y are erroneous
	 */
	public static double euclidian(double [] x, double [] y) throws IOException {
		return minkowski(2, x, y);
	}
	/**
	 * Finds the city block distance between x and y
	 * @param x : {double[]} the first element 
	 * @param y : {double[]} the second element 
	 * @return {double}
	 * @throws IOException if lenghts of x and y are erroneous
	 */
	public static double cityBlock(double [] x, double [] y) throws IOException {
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
	public static double minkowski(int power, double [] x, double [] y)throws IOException {
		//error checking
		errorCheck(x, y);
		
		double result = 0.0;
		for (int i=0; i < x.length; i++) {
			result += Math.pow(Math.abs(x[i] - y[i]), power);
		}
		return Math.pow(result, 1/power);
	}
	/**
	 * Finds the sup distance (max distance between features)
	 * @param x : {double[]} the first element 
	 * @param y : {double[]} the second element 
	 * @return {double}
	 * @throws IOException if lenghts of x and y are erroneous
	 */
	public static double sup(double [] x, double [] y) throws IOException {
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
}
