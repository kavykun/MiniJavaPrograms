import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class printPermutations implements Runnable {

	double countPermutations;
	double stddev = 0.0;
	double totalDistances = 0.0;
	double countDistances = 0.0;
	double meanDistances = 0.0;
	int[] str;
	Points2DLabel[] arrayPoints;
	double exhaustTime;

	public printPermutations(int[] s, Points2DLabel[] a) {

		this.str = s;
		this.arrayPoints = a;

	}

	@Override
	public void run() {
		// TODO Auto-generated method stub

		try {

			long start = System.nanoTime();
			printCombinations(this.str, 1, this.arrayPoints.length);
			printCost(this.str);
			compute();
			long end = System.nanoTime();

			exhaustTime = end - start;

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.out.println("The exhastive search time: " + exhaustTime
				+ " nanoseconds");
		System.out.println("The total distances: " + totalDistances);
		System.out.println("The mean distances: " + meanDistances);
		System.out.println("The standard devition: " + stddev);

	}

	/**************************************************************************
	 * Exhaustive Search Running time: O(N!)
	 **************************************************************************/

	/**
	 * Exhaust Search Second Method for permutation
	 * 
	 * @param str
	 *            the inputted string
	 * @param k
	 *            the start index
	 * @param n
	 *            the next index
	 * @throws IOException
	 */

	public void printCombinations(int[] str, int i, int n) throws IOException {

		PrintWriter out = null;
		int j;
		if (i == n) {
			out = new PrintWriter(new FileWriter("exhaustPermutations.txt",
					true));
			out.write(Arrays.toString(str));
			// printCost(str);
			countPermutations++;
			System.out.println(countPermutations);
			out.println();
			out.close();
		}
		for (j = i; j < n; j++) {

			swap(str, i, j);
			printCombinations(str, i + 1, n);
			swap(str, i, j); // backtrack

		}

	}// end print combinations

	/* Function to swap values at two pointers */
	public void swap(int[] str, int x, int y) {

		int temp;
		temp = str[x];
		str[x] = str[y];
		str[y] = temp;
	}

	/**************************************************************************
	 * Write the costs of each permutation to a text file
	 **************************************************************************/

	/**
	 * Takes the file that contain the permutation Computers each cost of travel
	 * for each permutation
	 * 
	 * @throws IOException
	 */
	public void printCost(int[] s) throws IOException {
		PrintWriter out = null;
		double totalDistanceTraveled = 0.0;
		double distance = 0.0;

		System.out.println("Printing..");

		Points2DLabel currentPoint2D = null;
		Points2DLabel previousPoint2D = null;

		previousPoint2D = arrayPoints[s[0]];

		for (int i = 1; i < s.length; i++) {

			currentPoint2D = arrayPoints[s[i]];
			distance = previousPoint2D.getPoint().distance(
					currentPoint2D.getPoint());

			previousPoint2D = currentPoint2D;

			totalDistanceTraveled += distance;

		}// end i for

		distance = previousPoint2D.getPoint().distance(
				arrayPoints[0].getPoint());

		totalDistanceTraveled += distance;

		out = new PrintWriter(new FileWriter("distanceCost.txt", true));
		out.println(totalDistanceTraveled + ",");
		out.close();
		totalDistanceTraveled = 0.0;
	}// end printCost

	/**************************************************************************
	 * Compute the distances of each permutation
	 **************************************************************************/

	public void compute() throws IOException {

		File file = null;
		double xisq = 0.0;
		double xisqSum = 0.0;
		double xiSumsq = 0.0;
		stddev = 0.0;

		file = new File("distanceCost.txt");

		long start2 = System.nanoTime();

		BufferedReader br = new BufferedReader(new FileReader(file));

		for (String line; (line = br.readLine()) != null;) {

			line = line.substring(0, line.length() - 1);

			double currentDistance = Double.parseDouble(line);

			totalDistances += currentDistance;
			countDistances++;

			xisq = Math.pow(currentDistance, 2);
			xisqSum += xisq;
			xiSumsq += currentDistance;

		}// end for
		br.close();
		long time2 = System.nanoTime() - start2;

		System.out.printf(
				"Took %.3f seconds to read to a %d MB file, rate: %.1f MB/s%n",
				time2 / 1e9, file.length() >> 20, file.length() * 1000.0
						/ time2);

		meanDistances = totalDistances / countDistances;

		// compute standard dev.

		xiSumsq = Math.pow(xiSumsq, 2);

		stddev = Math.sqrt((xisqSum - (xiSumsq) / countDistances)
				/ (countDistances - 1));

	}// end getValues

}
