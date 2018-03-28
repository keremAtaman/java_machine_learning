package genetic_algorithms;

import java.util.Random;

public class Examples {
	public class maze{
		public void maze() {
			int[][] maze = createMaze();
			Object[] possibleValues = {false, true};
			int numChromosomes = 64;
			Object[] temp_chromosome = Chromosomes.chromosomeCreator(numChromosomes, possibleValues);
			boolean[] chromosome = new boolean[numChromosomes];
			for (int i = 0; i < numChromosomes; i++) {
				chromosome[i] = (boolean) temp_chromosome[i];
			}
			
			
		}
		
		private int [][] createMaze(){
			int[] maze_temp = 
				   {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
					1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1,
					8, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1,
					1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
					1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
					1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
					1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1,
					1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 5,
					1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
					1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
			int[][] maze = new int [10][15];
			for (int i = 0; i < 10; i++) {
				for (int j = 0; j < 15; j++) {
					maze[i][j] = maze_temp[15 * i + j];
				}
			}
			return maze;
		}
		private int[] move(int[][] maze, int[] currentPosition, boolean[] directionGene) {
			int directionInt = Examples.bool2Int(directionGene);
			int[] direction = new int[2];
			switch (directionInt) {
				case 0:
					direction[0] = 1;
					direction[1] = 0;
					return moveHelper(maze, currentPosition, direction);
				case 1:
					direction[0] = -1;
					direction[1] = 0;
					return moveHelper(maze, currentPosition, direction);
				case 2:
					direction[0] = 0;
					direction[1] = -1;
					return moveHelper(maze, currentPosition, direction);
				case 3:
					direction[0] = 0;
					direction[1] = -1;
					return moveHelper(maze, currentPosition, direction);
				default:
					throw new IllegalArgumentException("directionGene must have 2 elements");
			}
			
		}
		private int[] moveHelper(int[][] maze, int[] currentPosition, int[] direction) {
			int[] position = new int[2];
			for (int i = 0; i < 2; i++) {
				position[i] = currentPosition[i] + direction[i];
			}
			
			if (moveHelperHelper(maze, position) == true) {
				return position;
			}else {
				return currentPosition;
			}
		}
		private boolean moveHelperHelper(int[][] maze, int[] position) {
			//ensure position does not go out of the board
			if ((position[0] > maze.length) | (position[1] > maze[0].length) | 
					(position[0] < 0) | (position[1] < 0)) {
				return false;
			}
			if (maze[position[0]][position[1]] == 1) {
				return false;
			}else{
				return true;
			}
		}
	}
	
	private static int bool2Int(boolean[] boolArray) {
		int result = 0;
		for (int i = 0; i < 4; ++i) {
	    	if (boolArray[i]) result |= (1 << i);
		}
		return result;
	}
}
