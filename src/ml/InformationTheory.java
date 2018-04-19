package ml;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class InformationTheory {//TODO use pmf instead of element occurences?
	public static double entropy(double[] pmf) {
		double entropy = 0.0;
		for (int i = 0; i <  pmf.length; i++) {
			if (pmf[i] != 0) {
				entropy -= pmf[i] * Math.log(pmf[i]);
			}
		}
		return entropy;
	}
	/**
	 * This is for 2 variables (namely x and y for explanation) only
	 * @param jointPmf row are x values, cols are y values
	 * @return
	 */
	public static double jointEntropy(double[][] jointPmf) {
		double jointEntropy = 0.0;
		for (int i = 0; i < jointPmf.length; i++) {
			for (int j = 0; j < jointPmf[0].length; j++) {
				if(jointPmf[i][j] != 0) {
					jointEntropy -= jointPmf[i][j] * Math.log(jointPmf[i][j]);
				}
			}
		}
		return jointEntropy;
	}
	public static double jointEntropy(Map<int[], Double> jointPmf) {
		double jointEntropy = 0.0;
		for (Entry<int[], Double> entry: jointPmf.entrySet()) {
			jointEntropy -= entry.getValue() * Math.log(entry.getValue());
		}
		return jointEntropy;
	}
	
	public static double mutualInformation(double[] x_alone, double[] y, double[][] x) {
		return entropy(x_alone) + entropy(y) - jointEntropy(x);
	}
	
	public static double mutualInformation(double[] x, double[] y, Map<int[], Double> jointPmf) {
		return entropy(x) + entropy(y) - jointEntropy(jointPmf);
	}
	
	public static double[][] featureClusters(double[] feature, int maxIterations, Distance.distanceMetric distanceMetric){
		double[][] x = new double[feature.length][1];
		for (int i = 0; i < feature.length; i++) {
			x[i][0] = feature[i];
		}
		return Clustering.ClusterEvaluation.optimalCenters(
				2, x.length, maxIterations, x, distanceMetric);
	}
	
	public static double[] featurePmf(double[][] featureCenters, double[] feature, Distance.distanceMetric distanceMetric){
		double[][] x = new double[feature.length][1];
		double[] pmf = new double[featureCenters.length];
		for (int i = 0; i < feature.length; i++) {
			x[i][0] = feature[i];
		}
		int[] memberships = Clustering.clusterMembership(featureCenters, x, distanceMetric);
		for (int i = 0; i < memberships.length; i++) {
			pmf[memberships[i]] += 1; 
		}
		for (int i = 0; i < memberships.length; i++) {
			pmf[i] = (double) pmf[i] / (double)memberships.length ;
		}
		
		return pmf;
	}
	
	public static double[] clusterPmf(int[] memberships, int numClusters) {
		double[] clusterPmf = new double[numClusters];
		for (int i = 0; i < memberships.length; i++) {
			clusterPmf[memberships[i]] += 1.0 / (double)memberships.length;
		}
		return clusterPmf;
	}
	
	
	private static boolean[] int2Bool(int b, int size) {
	    boolean [] result = new boolean[size];
		for (int i = 0; i < size; ++i) {
			result[i] = (b & (1 << i)) != 0;
		}
		return result;
	}
	
	private static double[] featureExtractor(double[][] x, int numFeature) {
		double[] result = new double[x.length];
		for (int i=0; i < x.length; i++) {
			result[i] = x[i][numFeature];
		}
		return result;
	}
	
	private static double[] jointPmf(int[] memberships, int[] featureMemberships) {
		double[][] placeholder = new double[memberships.length][featureMemberships.length];
		int counter = 0;
		for (int i = 0; i < memberships.length; i++) {
			for (int j = 0; j < featureMemberships.length; j++) {
				placeholder[memberships[i]][featureMemberships[j]]+= 1;
				counter ++;
			}
		}
		double[] result = new double[memberships.length * featureMemberships.length];
		for (int i = 0; i < memberships.length; i++) {
			for (int j = 0; j < featureMemberships.length; j++) {
				result[featureMemberships.length * i + j] = (double) (placeholder[i][j] / counter);
			}
		}
		return result;
	}
	
	private static double mRMRHelper(double[][] x, double[][] centers, int maxIterations, Distance.distanceMetric distanceMetric){
		//gets the stuff inside max bracket
		
		int[] memberships = Clustering.clusterMembership(centers, x, distanceMetric);
		//Calculate Entropy(Clusters)
		double clusterEntropy = InformationTheory.entropy(
				InformationTheory.clusterPmf(memberships, centers.length));
		
		//To calculate Entropy(Features), find the feature pmfs. 
		
		//featurePmf[i] is the ith feature's pmf
		List<double[]> featurePmfList = new ArrayList<double[]>();
		List<double[][]>featureCentersList = new ArrayList<double[][]>();
		double[] featureEntropyArray = new double[x[0].length];
		double[][] featureExtract = new double[x.length][1];
		double [] placeholder = new double[x.length];
		for (int i = 0; i < x[0].length; i++) {
			placeholder = InformationTheory.featureExtractor(x, i);
			for (int k = 0; k < placeholder.length; k++) {
				featureExtract[k][0] = placeholder[k];
			}
			//find the centers for feature i
			featureCentersList.add(Clustering.ClusterEvaluation.optimalCenters(
					1, x.length, maxIterations, featureExtract, distanceMetric)) ;
			//add the ith feature's pmf
			featurePmfList.add(
					InformationTheory.featurePmf(
							featureCentersList.get(i), featureExtract[0], distanceMetric));
			//Get the entropy for the ith feature
			featureEntropyArray[i] = InformationTheory.entropy(featurePmfList.get(i));
		}
		
		
		double sum = 0.0;
		double sum2 = 0.0;
		double jointEntropy = 0.0;	
		for (int i=0; i < x[0].length; i++) {
			//Get joint entropy of f and c
			
			//get the feature
			placeholder = InformationTheory.featureExtractor(x, i);
			for (int k = 0; k < placeholder.length; k++) {
				featureExtract[k][0] = placeholder[k];
			}
			int[] featureMemberships = Clustering.clusterMembership(
					featureCentersList.get(i), featureExtract, distanceMetric);
			jointEntropy = InformationTheory.entropy(InformationTheory.jointPmf(
					memberships, featureMemberships));
			//Get the MI(f,c)
			sum += clusterEntropy + featureEntropyArray[i] - jointEntropy;
			
			for (int j = 0; j < x[0].length; j++) {
				placeholder = InformationTheory.featureExtractor(x, j);
				for (int k = 0; k < placeholder.length; k++) {
					featureExtract[k][0] = placeholder[k];
				}
				int [] featureMemberships2 = Clustering.clusterMembership(
						featureCentersList.get(j), featureExtract, distanceMetric);
				jointEntropy = InformationTheory.entropy(InformationTheory.jointPmf(
						featureMemberships2, featureMemberships));
				sum2+= featureEntropyArray[i] + featureEntropyArray[j] - jointEntropy;
			}
		}
		return sum / (double) x[0].length - sum2 / (double) Math.pow(x[0].length, 2);
	}

	public static boolean[] mRMR(double[][] x, double[][] centers, int maxIterations, Distance.distanceMetric distanceMetric){
		boolean[] featuresToRemove = int2Bool(0, x[0].length);
		boolean[] bestFeaturesToRemove = new boolean[x[0].length];
		int counter = 0;
		int counter2 = 0;
		int numOnes = 0;
		double max = 0.0;
		double score = 0.0;
		int numFeatures = x[0].length;
		
		//Go over every element deletion possibility
		for (counter = 0; counter < Math.pow(numFeatures, 2) - 1; counter++) {
			featuresToRemove = int2Bool(counter, x[0].length);
			numOnes = 0;
			//count the # of trues in features to remove
			for (int i = 0; i < featuresToRemove.length; i++) {
				if (featuresToRemove[i] == true) {
					numOnes ++;
				}
			}
			
			//create a new array with features to be removed being removed
			double[][] newX = new double[x.length][x[0].length - numOnes];
			double[][] newCenters = new double[centers.length][x[0].length - numOnes];
			counter2 = 0;
			for (int j = 0; j < x[0].length; j++) {
				//if this column is to be skipped, increase the count for the number of columns to skip
				if (featuresToRemove[j]==true) {
					counter2 ++;
				}
				//fill the column if it is not to be removed
				else {
					for (int i=0; i < x.length; i++) {
						//if the feature would not be removed, add it to 
						newX[i][j - counter2] = x[i][j];
						
					}
					for (int k = 0; k < centers.length; k++) {
						newCenters[k][j - counter2] = centers[k][j];
					}
				}
			}
			
			score = InformationTheory.mRMRHelper(newX, newCenters, maxIterations, distanceMetric);
			
			//if score is NaN
			if (Double.isNaN(score)) {
				throw new ArithmeticException("score is NaN");
			}
			//If the current score is greater than the max, current score and features removed
			//becomes the new best
			if(score > max) {
				bestFeaturesToRemove = featuresToRemove;
				max = score;
			}
		}
		return bestFeaturesToRemove;
	}
}
