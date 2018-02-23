package ml;

import java.io.IOException;
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
	public static double covariance(double[] x, double[] y) throws IOException {
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
		
		double meanX = mean(x);
		double meanY = mean(y);
		double result = 0.0;
		for (int i = 0; i < x.length; i++) {
			result += (x[i] - meanX) * (y[i] - meanY);
		}
		return result / x.length;
	}
	public static double[][] covarianceMatrix(double[][]x) throws IOException {
		double[][] result = new double[x.length][x[0].length];
		for (int i = 0; i < x.length; i++) {
			for (int j = 0; j < x[0].length; j++) {
				result[i][j] = Statistics.covariance(x[i], x[j]);
			}
		}
		return result;
	}
	public static double pearsonCorrelationCoefficient(double[] x, double [] y) throws IOException {
		return Statistics.covariance(x, y) / (Statistics.stdev(x) * Statistics.stdev(y));
	}
}
