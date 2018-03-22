package ml;

public class Dataset {
	public static class Generate{
		/**
		 * Creates numElements / elementsPerCluster clusters
		 * @param elementsPerCluster how many elements should there be in each cluster?
		 * @param numElements number of elements total in the dataset
		 * @param numFeatures number of features in dataset
		 * @param stdevFeatures standard deviation between features
		 * @param stdevClusters standard deviation between clusters
		 * @return
		 */
		public static double[][] clusters(int elementsPerCluster, int numElements, int numFeatures, double stdevFeatures, 
				double stdevClusters){
			double[][] result = new double[numElements][numFeatures];
			for (int i = 0; i < numElements; i++) {
				for (int j = 0; j < numFeatures; j++) {
					result[i][j] = Math.random() * stdevFeatures + stdevClusters * 
							Math.floor(i / elementsPerCluster);
				}
			}
			return result;
		}
	}
}
