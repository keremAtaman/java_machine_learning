package genetic_algorithms;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Random;

public class Chromosomes {
	
	/**
	 * Different methods for cross over
	 */
	public enum crossoverMethod{
		/**
		 * At the crossover point, if crossover rate check is passed, flip the chromosomes of parents
		 */
		CROSSOVER_POINT
	}
	public enum mutationMethod{
		/**
		 * Randomly mutates elements
		 */
		BIT_STRING,
		/**
		 * Mutates all elements
		 */
		FLIP_BIT
	}
	
	/**
	 * Crosses over x and y according to method specified
	 * @param x
	 * @param y
	 * @param crossoverRate
	 * @param crossoverMethod
	 * @return a 2 by x.length array that has the children in it's rows
	 * @throws IOException
	 */
	public static Object[][] crossover(Object[] x, Object [] y, double crossoverRate,
			Chromosomes.crossoverMethod crossoverMethod) throws IOException {
		if (crossoverRate >1 | crossoverRate < 0) {
			throw new IOException("crossover rate must be between 0 and 1, inclusive");
		}
		if (x == null){
			throw new IOException("x must be non-empty");
		}
		if (x.length != y.length) {
			throw new IOException ("Length of x and y must be the same");
		}
		
		
		switch (crossoverMethod) {
		case CROSSOVER_POINT:
			return noRestrictionCrossover(x, y, crossoverRate);
		default:
			throw new IOException("Invalid crossover method");
		}
	}
	
	private static Object[][] noRestrictionCrossover(Object[] x, Object [] y, double crossoverRate) throws IOException {
		
		Object[][] result = new Object[2][x.length];
		Random rng = new Random();
		// do a crossover
		if (rng.nextFloat() <= crossoverRate) {
			//determine the point of cross over
			int crossoverPoint = rng.nextInt(x.length);
			for (int i = 0; i < x.length; i++) {
				if (i < crossoverPoint) {
					result[0][i] = x[i];
					result[0][i] = y[i];
				}else {
					result[0][i] = y[i];
					result[0][i] = x[i];
				}
			}
		//no crossover
		}else {
			result[0] = x;
			result[1] = y;
		}
		return result;
	}
	/**
	 * 
	 * @param x
	 * @param mutationRate
	 * @param mutationMethod
	 * @param mutationDict has the format element to mutate -> possible mutations for the element
	 * @return
	 * @throws IOException
	 */
	public static Object[] mutate(Object[] x, double mutationRate, mutationMethod mutationMethod, Map <Object, ArrayList<Object>> mutationDict) throws IOException {
		if (mutationRate >1 | mutationRate < 0) {
			throw new IOException("mutation rate must be between 0 and 1, inclusive");
		}
		if (mutationDict == null) {
			throw new IOException("mutation dictionary must be filled");
		}
		if (x==null) {
			throw new IOException("x must be non-empty");
		}
		switch(mutationMethod) {
		case BIT_STRING:
			return bitStringMutation(x, mutationRate, mutationDict);
		case FLIP_BIT:
			return flipBitMutation(x, mutationRate, mutationDict);
		default:
			return null;
		}
		
	}
	
	private static Object[] bitStringMutation(Object[] x, double mutationRate, Map <Object, ArrayList<Object>> mutationDict) {
		Random rng = new Random();
		int dictSelector = 0;
		for (int i = 0; i < x.length; i++) {
			//mutate
			if (rng.nextFloat() < mutationRate) {
				dictSelector = mutationDict.get(x[i]).size();
				x[i] = mutationDict.get(x[i]).get(dictSelector);
			}
		}
		return x;
	}
	private static Object[] flipBitMutation(Object[] x, double mutationRate, Map <Object, ArrayList<Object>> mutationDict) {
		Random rng = new Random();
		int dictSelector = 0;
		//if mutation test is passed, mutate every element
		if (rng.nextFloat() < mutationRate) {
			for (int i = 0; i < x.length; i++) {
				dictSelector = mutationDict.get(x[i]).size();
				x[i] = mutationDict.get(x[i]).get(dictSelector);
			}
		}
		return x;
	}
}
