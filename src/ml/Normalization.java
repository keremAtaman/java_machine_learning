package ml;

public class Normalization {
	/**
	 * 
	 * @param x rows are entries, columns are features
	 * @return standard score normalized version of x
	 */
	public static double[][] zScore(double[][] x){
		if (x == null) {
			throw new IllegalArgumentException("x has to be nonempty");
		}
		
		double mean = 0.0;
		double stdev = 0.0;
		
		double[][] result = new double[x.length][x[0].length];
		//an array where the feature column will be stored
		double[] featureArray = new double[x[0].length];
		
		for (int j = 0; j < x[0].length; j++) {
			//fill the feature array
			for (int i = 0; i < x.length; i++) {
				featureArray[i] = x[i][j];
			}
			//find mean and stdev of the feature
			mean = Statistics.mean(featureArray);
			stdev = Statistics.stdev(featureArray);
			//prevent 0 cases from causing trouble
			if (stdev == 0) {
				stdev = 1;
			}
			
			for (int i = 0; i < x.length; i++) {
				result[i][j] = (featureArray[i] - mean) / stdev;
			}
		}
		
		return result;
	}
}
