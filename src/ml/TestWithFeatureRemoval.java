package ml;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TestWithFeatureRemoval {

	public static void main(String[] args) throws IOException {
		//TODO: Fix all features being removed AND identify the names of features being removed
		//TODO: Fix x is empty bug
		//TODO: +-1 k using last center values to start
		String location = 
				"C:/Users/K.Ataman/Dropbox (TTS)/TTS Development/Machine Learning/Cpty POC/QT Mar 20";
		String[] header = {"NUM_CCY_TRADED","AVG_TRADE_SIZE" ,"AVG_TRADE_VOLUME_MNTHLY",
				 "CONSISTENCY_INDEX","%_TRADES_BEFORE_NEWS","%_TRADES_AFTER_NEWS" ,
				"%_TRADES_DURING_HIGH_VOLATILITY" ,"FILL_RATIO" ,"%_ROUND_NUMBER_TRADED" ,"%_SPOT", "%_PERIODICITY", "%_USDCAD"};	
		//reads the csv and populates x accordingly
		String[][] stringX =  CSVUtilities.csv_reader_with_header(location + "/aggregationsFinal.csv", header.length, 39, header);
		double[][]x = new double[stringX.length][stringX[0].length];
		for (int i = 0; i < stringX.length; i++) {
			for (int j = 0; j < stringX[0].length; j++) {
				x[i][j] = Double.parseDouble(stringX[i][j]);
			}
		}
		
		Distance.distanceMetric distanceMetric = Distance.distanceMetric.EUCLIDIAN;
		
		double score = 0.0;
		//ensure we get a good score
		double minScore = 0.5;
		int maxIterations = 300;
		List<double[]> centersHolder = new ArrayList<double []>();
		while(score < minScore) {
			//get the results given the current features
			double[][][] result = Clustering.ClusterEvaluation.clusteringProcess(x, maxIterations, distanceMetric);
			
			//get the centers and store them in a double[][]
			centersHolder.clear();
			for (int i = 0; i < result[0].length; i++) {
				if (result[0][i]!= null) {
					centersHolder.add(result[0][i]);
				}
			}
			double [][] centers = new double[centersHolder.size()][centersHolder.get(0).length];
			for (int i = 0; i < centersHolder.size(); i++) {
				for (int j = 0; j < centersHolder.get(0).length; j++) {
					centers[i][j] = centersHolder.get(i)[j];
				}
			}
			
			//calculate result details
			double[][] newX = result[1];
			int[] clusterMembership = Clustering.clusterMembership(centers, newX, distanceMetric);
			
			score = Clustering.ClusterEvaluation.silhouetteMethod(centers, newX, distanceMetric);
			//we have achieved a good score, do the bookkeeping
			if (score > minScore) {
				System.out.println("The optimal number of clusters is:");
				System.out.println(centers.length);
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
	}

}
