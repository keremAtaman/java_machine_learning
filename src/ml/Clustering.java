package ml;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class Clustering {
	
	/**
	 * Calculates distance between x and y using the given distance metric
	 * @param x
	 * @param y
	 * @param distanceMetric
	 * @return
	 * @throws IOException
	 */
	private static double calculateDistance(double[] x, double[] y, Distance.distanceMetric distanceMetric
			) throws IOException {
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
	/**
	 * calculates which cluster each element of x belongs to
	 * @param centers {double[numClusters][numFeatures]} a matrix that stores the centers of each cluster
	 * @param x {double[numElements][numFeatures]} the data that got clustered
	 * @param distanceMetric {Distance.distanceMetric} 
	 * @return {int[x.length]} An array that shows which cluster each element of x belongs to. ie, 
	 * output[0] = 2 means the first element of x belongs to the third cluster
	 * @throws IOException
	 */
	public static int [] clusterMembership(double[][] centers, 
			double x[][], Distance.distanceMetric distanceMetric) throws IOException{
		int [] clusterMembers = new int[x.length];
		double distance;
		double minDistance = 0.0;
		int closestCluster = -1;
		// calculate which cluster each element belongs to
		for (int i = 0; i < x.length ; i++) {
			for (int j = 0; j < centers.length; j++) {
				distance = Clustering.calculateDistance(x[i], centers[j], distanceMetric);
				if (j==0) {
					minDistance = distance;
					closestCluster = 0;
				}else if (distance < minDistance){
					minDistance = distance;
					closestCluster = j;
				}
			}
			clusterMembers[i] = closestCluster;
		}
		return clusterMembers;
	}
	public static class ClusterEvaluation{
		/**
		 * Almost there, but for some reason rearrangedMembers is not filled properly
		 * @param centers
		 * @param x
		 * @param distanceMetric
		 * @return
		 * @throws IOException
		 */
		public static double silhouetteMethod(double[][] centers, double[][] x, Distance.distanceMetric distanceMetric) throws IOException {
			double score = 0.0;
			//Get which cluster each element belongs to
			int[] clusterMemberships = Clustering.clusterMembership(centers, x, distanceMetric);
			//first index is clusters, second index is elements
			List<List<double[]>> rearrangedMembers= new ArrayList<List<double[]>>();
			List<double[]> clusterMembers = new ArrayList<double[]>();
			
			//populate the rearrangedMembers List
			for(int i=0; i < centers.length; i++) {
				for (int j = 0; j < x.length; j++) {
					if(clusterMemberships[j] == i) {
						clusterMembers.add(x[j]);
						//System.out.println(clusterMembers.get(0)[0]);
					}
				}
				rearrangedMembers.add(clusterMembers);
				//System.out.println(rearrangedMembers.get(0).get(0)[0]);
				clusterMembers.clear();
			}
			System.out.println(rearrangedMembers.size());
			System.out.println(rearrangedMembers.get(0).size());
			//calculate within and inter cluster distances for all elements
			//within cluster score
			double SSW = 0.0;
			//inter cluster score
			double SSB = 0.0;
			for (int i = 0; i < rearrangedMembers.size(); i++) {
				for (int j = 0; j < rearrangedMembers.get(i).size(); j++) {
					SSW = withinClusterDistance(j, rearrangedMembers.get(i), distanceMetric);
					SSB = Clustering.ClusterEvaluation.interClusterDistance(i, rearrangedMembers, distanceMetric);
					score += (SSW - SSB) / Math.max(SSB, SSW);
					System.out.println(SSW);
					System.out.println(SSB);
				}
			}
			return score / (double)x.length;
		}
		
		private static double withinClusterDistance(int elementNumber, List<double[]> clusterMembers, Distance.distanceMetric distanceMetric) throws IOException {
			double score = 0.0;
			for(int i = 0; i < clusterMembers.size(); i++) {
				if(i != elementNumber) {
					score += Clustering.calculateDistance(clusterMembers.get(elementNumber), clusterMembers.get(i), distanceMetric);
				}
			}
			return score / (double)clusterMembers.size();
		}
		
		private static double interClusterDistance(int clusterNumber, List<List<double[]>> rearrangedMembers, Distance.distanceMetric distanceMetric) throws IOException {
			double score = 0.0;
			for (int i = 0; i < rearrangedMembers.size(); i++) {
				if (i!= clusterNumber) {
					score+= Clustering.ClusterEvaluation.withinClusterDistance(-1, rearrangedMembers.get(i), distanceMetric);
				}
			}
			return score / (double)rearrangedMembers.size();
		}
		
		public static int elbowMethod(double[] listOfScores) {
			double[] secondDerivativeList = new double[listOfScores.length];
			double maxScore = 0.0;
			int maxIndex = 0;
			for (int i = 0; i < listOfScores.length - 1; i++) {
				secondDerivativeList[i] = Math.abs(listOfScores[i+1] + listOfScores[i-1] - 2 * listOfScores[i]);
				if (secondDerivativeList[i] > maxScore) {
					maxIndex = i;
				}
			}
			return maxIndex;
		}
	}
	
	
	
	public static int [] clusterMembership(List<double[]> centers, 
			double x[][], Distance.distanceMetric distanceMetric) throws IOException{
		int [] clusterMembers = new int[x.length];
		double distance;
		double minDistance = 0.0;
		int closestCluster = -1;
		// calculate which cluster each element belongs to
		for (int i = 0; i < x.length ; i++) {
			for (int j = 0; j < centers.size(); j++) {
				distance = Clustering.calculateDistance(x[i], centers.get(j), distanceMetric);
				if (j==0) {
					minDistance = distance;
					closestCluster = 0;
				}else if (distance < minDistance){
					minDistance = distance;
					closestCluster = j;
				}
			}
			clusterMembers[i] = closestCluster;
		}
		return clusterMembers;
	}
	
	private static double[][] initialCenters(int k, double [][] x){
		//get the min and max of each feature
		double [] min = new double[x[0].length];
		double [] max = new double[x[0].length];
		for (int i = 0; i < x.length; i++) {
			for (int j = 0; j < x[0].length; j++) {
				if (i==0) {
					min[j] = x[i][j];
					max[j] = x[i][j];
				}else if (min[j] > x[i][j]) {
					min[j] = x[i][j];
				}else if (max[j] < x[i][j]) {
					max[j] = x[i][j];
				}
			}
		}
		//holds the centers
		double[][] centers = new double[k][x[0].length];
		//holds the centers of last step for comparison
		//generate initial centers, make sure they don't exceed maximum and minimum of their
		//respective features
		for (int i=0; i < k; i++) {
			for (int j = 0; j < x[0].length; j++) {
				centers[i][j] = Math.random() * (max[j] - min[j]) + min[j];
			}
		}
		return centers;
	}
	
	
	/**
	 * Does knn on a given dataset
	 * @param k number of clusters
	 * @param x data to be clustered, where rows are entries and columns are features
	 * @param initialCenters if not empty, gives the initial centers to be used
	 * @param distanceMetric distance metric to be used when calculating distance
	 * @return the centers of each cluster
	 * @throws IOException
	 */
	public static double[][] kMeans(int k, double[][] x, double[][] initialCenters, Distance.distanceMetric distanceMetric
			) throws IOException {
		//error checking
		if (k < 2) {
			throw new IOException("k has to be an intager greter than 1");
		}
		if (x.length < 2) {
			throw new IOException("x has to have at least two entries");
		}
		if (k >= x.length) {
			throw new IOException("k has to be less than the number of entries in x");
		}
		
		//holds the centers. Give initial random values to the center coordinates if no default value is given
		double[][] centers = initialCenters;
		if (initialCenters == null) {
			centers = Clustering.initialCenters(k, x);
		}
		//holds the centers of last step for comparison
		double[][] previousCenters = new double[k][x[0].length];
		
		boolean terminationCondition = false;
		//Holds who belongs to which cluster
		int []clusterMembers = new int[x.length];
		//holds how many people are there in each cluster
		int[] numClusterMembers = new int[k];
		int closestCluster = -1;
		double minDistance = 0.0;
		double distance = 0.0;
		
		while(terminationCondition == false) {
			previousCenters = centers;
			//reset groupMemberCount
			for (int i=0; i < x.length; i++) {
				clusterMembers[i] = 0;
			}
			for (int i=0; i < k; i++) {
				
			}
			
			// calculate which cluster each element belongs to
			for (int i = 0; i < x.length ; i++) {
				for (int j = 0; j < k; j++) {
					distance = Clustering.calculateDistance(x[i], centers[j], distanceMetric);
					if (j==0) {
						minDistance = distance;
						closestCluster = 0;
					}else if (distance < minDistance){
						minDistance = distance;
						closestCluster = j;
					}
				}
				clusterMembers[i] = closestCluster;
				numClusterMembers[closestCluster] += 1;
			}
			// reset centers
			for (int i=0; i < k; i++) {
				for (int j = 0; j < x[0].length; j++) {
					centers[i][j] = 0;
				}
			}
			//calculate new centers
			for (int j = 0; j < k; j++) {
				for (int l = 0; l < x.length; l++) {
					if (clusterMembers[l] == j) {
						for (int feature = 0; feature < x[0].length; feature++) {
							centers[j][feature] += x[l][feature]/numClusterMembers[j];
						}
					}
				}
			}
			//compare new centers with the old ones. If they are the same, break
			if (centers == previousCenters) {
				break;
			}
		}
		return centers;
	}
	/**
	 * INCOMPLETE: NEED TO DO MERGING STEP AND CLEAR OUT TO ENSURE ONLY List <double[]> centers is used
	 * @param k
	 * @param minElements
	 * @param maxNumberOfIterations
	 * @param maxVariance
	 * @param minDistance
	 * @param x
	 * @param distanceMetric
	 * @return
	 * @throws IOException
	 */
	public static double[][] ISODATA(int k, int minElements, int maxNumberOfIterations,
			double maxVariance, double minDistance,
			double[][] x, Distance.distanceMetric distanceMetric) throws IOException{
		
		//error checking
		if (k < 2) {
			throw new IOException("k has to be an intager greter than 1");
		}
		if (x.length < 2) {
			throw new IOException("x has to have at least two entries");
		}
		if (k >= x.length) {
			throw new IOException("k has to be less than the number of entries in x");
		}
		
		
		//holds the centers
		List<double[]> centers = new ArrayList<double[]>();
		//set initial centers
		double[][] temp = Clustering.initialCenters(k, x);
		for (int i = 0; i < temp.length; i++) {
			centers.add(temp[i]);
		}
		for(int iteration = 0; iteration < maxNumberOfIterations; iteration++) {
			//see where each element is in centers
			int[] clusterMembers = Clustering.clusterMembership(centers, x, distanceMetric);
			int[] clusterPopulation = new int[clusterMembers.length];
			//recieve the number of members for each cluster
			for (int i =0; i < clusterMembers.length; i++) {
				clusterPopulation[clusterMembers[i]]+=1;
			}
			
			//remove clusters with too few elements in it
			int i = 0;
			Iterator<double[]> it = centers.iterator();
			while (i < centers.size()) {
				
				double[] current = it.next();
				if (clusterPopulation[i] < minElements) {
					it.remove();
				}else {
					i++;
				}
			}
			double[][]newCenters = new double[centers.size()][x[0].length];
			Iterator<double[]> itrtr = centers.iterator();
			for (int a = 0; i < centers.size(); i++) {
				for (int b = 0; b < x[0].length; b++) {
					newCenters[a][b] = itrtr.next()[b];
				}
			}
			clusterMembers = Clustering.clusterMembership(newCenters, x, distanceMetric);
			int[] newClusterPopulation = new int[clusterMembers.length];
			
			
			//recieve the number of members for each cluster
			for (i =0; i < clusterMembers.length; i++) {
				newClusterPopulation[clusterMembers[i]]+=1;
			}
			//assign new centers
			for (i=0; i < clusterMembers.length; i++) {
				for (int j=0; j < x[0].length; j++) {
					newCenters[clusterMembers[i]][j] = x[i][j] / newClusterPopulation[clusterMembers[i]];
				}
			}
			//create covariance matrix
			double[][] covarianceMatrix = Statistics.covarianceMatrix(newCenters);
			
			
			//too few elements
			if (newCenters.length <= k / 2.0) {
				//find max covaraince of each cluster
				double[] maxCovariance = new double[covarianceMatrix.length];
				for (i=0; i < covarianceMatrix.length; i++) {
					for (int j = 0; j < covarianceMatrix[0].length; j++) {
						if (j==0) {
							maxCovariance[i] = covarianceMatrix[i][j];
						}else if (covarianceMatrix[i][j] > maxCovariance[i]) {
							maxCovariance[i] = covarianceMatrix[i][j];
						}
					}
				}
				
				//check if the variance exceeds the desired amount of variance
				Iterator<double[]> it2 = centers.iterator();
				double[] temp1 = new double[newCenters[0].length];
				double[] temp2 = new double[newCenters[0].length];
				for(i=0; i < covarianceMatrix.length; i++) {
					double[] current = it2.next();
					//there is enough variance and elements to split the cluster into two
					if ((maxCovariance[i] * maxCovariance[i] > maxVariance * maxVariance
							) & (newClusterPopulation[i] > 2 * minElements)) {
						for (int j = 0; j < covarianceMatrix[0].length; j++) {
							temp1[j] = current[j] + maxCovariance[j];
							temp2[j] = current[j] - maxCovariance[j];
						centers.add(temp1);
						centers.add(temp2);
						//finally, remove the splitted cluster
						it2.remove();
						}
					}
				}
				
				
			//too many elements	
			}else if (centers.size() >= 2*k) {
				//merge, see http://fourier.eng.hmc.edu/e161/lectures/classification/node12.html
			}
			
			
		}
		double[][]result = new double[centers.size()][x[0].length];
		for (int i = 0; i < result.length; i++) {
			for (int j = 0; j < result[0].length; j++) {
				result[i][j] = centers.get(i)[j];
			}
		}
		return result;
	}
}
