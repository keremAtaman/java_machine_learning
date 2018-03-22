package ml;

import java.io.IOException;

public class Test {
	public static final void main(String[] args) throws IOException {
		double [][] x = {{0.1, 0.2}, {0.11, 0.21}, {4, 4.1}, {4.5, 5.1}};
		Distance.distanceMetric distanceMetric = Distance.distanceMetric.EUCLIDIAN;
		double [][] centers = Clustering.ClusterEvaluation.optimalCenters(
				2, x.length, x, distanceMetric);
		int[] clusterMembership = Clustering.clusterMembership(centers, x, distanceMetric);
		boolean[] featuresToRemove = InformationTheory.mRMR(x, centers, distanceMetric);
		System.out.println("The optimal number of clusters is:");
		System.out.println(centers.length);
		System.out.println("here are cluster memberships:");
		for (int i=0; i < x.length; i++) {
			System.out.println(clusterMembership[i]);
		}
		System.out.println("Here is the Silhouette score:");
		System.out.println(Clustering.ClusterEvaluation.silhouetteMethod(centers, x, distanceMetric));
	}
}
