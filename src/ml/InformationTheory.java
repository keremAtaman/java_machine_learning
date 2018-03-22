package ml;

public class InformationTheory {
	public static double entropy(double[] elementOccurences) {
		int sum = 0;
		double entropy = 0.0;
		for (int i = 0; i < elementOccurences.length; i++) {
			sum += elementOccurences[i];
		}
		double pmf = 0;
		for (int i = 0; i <  elementOccurences.length; i++) {
			pmf = elementOccurences[i] /sum;
			entropy -= pmf * Math.log(pmf);
		}
		return entropy;
	}
	
	public static double jointEntropy() {
		
	}
}
