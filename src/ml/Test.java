package ml;

import java.io.IOException;

public class Test {
	public static final void main(String[] args) throws IOException{
		String location = 
				"C:/Users/K.Ataman/Dropbox (TTS)/TTS Development/Machine Learning/Cpty POC/QT Mar 20";
		String[] header = {"NUM_CCY_TRADED" ,"NUM_TRADES" ,"AVG_TRADE_SIZE" ,"AVG_TRADE_VOLUME_MNTHLY"
				,"AVG_TRADE_COUNT_MNTHLY" ,"CONSISTENCY_INDEX","%_TRADES_BEFORE_NEWS","%_TRADES_AFTER_NEWS" ,
				"%_TRADES_DURING_HIGH_VOLATILITY" ,"FILL_RATIO" ,"%_ROUND_NUMBER_TRADED" ,"%_SPOT","%_FORWARDS", "%_PERIODICITY"};
		String[][] stringX =  CSVUtilities.csv_reader_with_header(location + "/aggregationsFinal.csv", header.length, 39, header);
		double[][]x = new double[stringX.length][stringX[0].length];
		for (int i = 0; i < stringX.length; i++) {
			for (int j = 0; j < stringX[0].length; j++) {
				x[i][j] = Double.parseDouble(stringX[i][j]);
			}
		}
		//double [][] x = {{0.1, 0.2}, {0.11, 0.21}, {4, 4.1}, {4.5, 5.1}};
		Distance.distanceMetric distanceMetric = Distance.distanceMetric.EUCLIDIAN;
		
		
		double [][] centers = Clustering.ClusterEvaluation.optimalCenters(
				2, x.length, x, distanceMetric);
		int[] clusterMembership = Clustering.clusterMembership(centers, x, distanceMetric);
		boolean[] featuresToRemove = InformationTheory.mRMR(x, centers, distanceMetric);
		double[][] newX = x;
		
		/**
		double[][][] results = Clustering.ClusterEvaluation.clusteringProcess(x, distanceMetric);
		double[][] centers = results[0];
		double[][] newX = results[1];
		int[] clusterMembership = Clustering.clusterMembership(centers, newX, distanceMetric);
		**/
		
		
		System.out.println("The optimal number of clusters is:");
		System.out.println(centers.length);
		/**
		 * Outputs cluster memberships. Use for small sets only
		System.out.println("here are cluster memberships:");
		for (int i=0; i < newX.length; i++) {
			System.out.println(clusterMembership[i]);
		}
		**/
		System.out.println("Here is the Silhouette score:");
		System.out.println(Clustering.ClusterEvaluation.silhouetteMethod(centers, newX, distanceMetric));
		System.out.println("Here is the amount of features removed:");
		System.out.println(x[0].length - newX[0].length);
		String [] clusterMembershipString= new String[clusterMembership.length];
		for (int i = 0; i < clusterMembership.length; i++) {
			clusterMembershipString[i] = String.valueOf(clusterMembership[i]);
		}
		String[] derp = {"cluster_id"};
		CSVUtilities.csv_writer_with_header(location, "/clusters.csv",derp ,clusterMembershipString);
		
		String[] headers = new String [centers[0].length];
		for (int i = 0; i < headers.length; i++) {
			headers[i] = Integer.toString(i);
		}
		CSVUtilities.csv_writer_with_header(location, "/cluster_details.csv",headers ,centers);
	}
}
