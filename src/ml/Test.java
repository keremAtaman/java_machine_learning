package ml;

public class Test {
	public static final void main(String[] args){
		double [][] x = {{0.1, 0.2}, {0.11, 0.21}, {4, 4.1}, {4.5, 5.1}};
		Distance.distanceMetric distanceMetric = Distance.distanceMetric.EUCLIDIAN;
		/**
		 * The below is for testing without feature removal
		double [][] centers = Clustering.ClusterEvaluation.optimalCenters(
				2, x.length, x, distanceMetric);
		int[] clusterMembership = Clustering.clusterMembership(centers, x, distanceMetric);
		boolean[] featuresToRemove = InformationTheory.mRMR(x, centers, distanceMetric);
		**/
		double[][][] results = Clustering.ClusterEvaluation.clusteringProcess(x, distanceMetric);
		double[][] centers = results[0];
		double[][] newX = results[1];
		int[] clusterMembership = Clustering.clusterMembership(centers, newX, distanceMetric);
		System.out.println("The optimal number of clusters is:");
		System.out.println(centers.length);
		System.out.println("here are cluster memberships:");
		for (int i=0; i < newX.length; i++) {
			System.out.println(clusterMembership[i]);
		}
		System.out.println("Here is the Silhouette score:");
		System.out.println(Clustering.ClusterEvaluation.silhouetteMethod(centers, newX, distanceMetric));
	}
}
