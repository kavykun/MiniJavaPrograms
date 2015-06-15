import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Test to manipulate Eigenvalues, Eigenvectors, and the Genetic Algorithm.
 * 
 * @author Kavy Rattana
 * @date 03/19/2014
 * 
 */
public class test2part2 {

	/**************************************************************************
	 * Variables
	 * 
	 **************************************************************************/

	public char[] labelsArray; // array to hold on the labels
	public String labelsString, fileName, permutationStrings;
	public LinkedList<String> permutations = new LinkedList<String>();
	public LinkedList<String> randomPermutations = new LinkedList<String>();
	public LinkedList<Permutation> permutationList = new LinkedList<Permutation>();
	public LinkedList<Permutation> childPopulation;
	public Points2DLabel[] arrayPoints;
	public int choice;
	public double[] distanceArray;
	public double meanDistances = 0.0, stddev = 0.0, exhaustTime = 0.0,
			randomSearchTime = 0.0, geneticAlgorithmSearchTime = 0.0;
	double totalDistances = 0.0;
	double countDistances = 0.0;
	double totalOverallDistances = 0.0;
	double countPermutations = 0.0;
	double totalfitness = 0.0;
	double xi = 0.0;
	double xisq = 0.0;
	double xisqSum = 0.0;
	double xiSumsq = 0.0;
	public Permutation[] mainPopulation; // the main population
	public int topIndividuals = 5;
	public int populationSize; // size of the population
	public int numberOfTests; // max number of iterations
	double mutationRate; // probability of mutation
	double crossoverRate; // probability of crossover
	public static double uniformRate = 0.5;
	private static Random m_rand = new Random(); // randomizer
	public double[] fitnesses; // array for top5 individual's fitness
	public Permutation[] top5; // array for top5 individuals
	public Permutation[] bestIndividuals;
	public int[] val;
	public int now = 0;
	public double distanceTotal = 0.0;
	public double sumXSquared;
	public double binWidth, max, min;
	public int V;
	public int numBins = 150;
	public double[] bins = new double[numBins];
	public long count = 0;
	public double countBinAdds;

	/**************************************************************************
	 * Running the program The text file should be called "data2.txt"
	 **************************************************************************/

	public static void main(String[] args) throws IOException {

		test2part2 obj = new test2part2();
		obj.run();

	}// end main

	public void run() throws IOException {

		createMatrix();
		double xi = 0.0;
		double xisq = 0.0;
		double xisqSum = 0.0;
		double xiSumsq = 0.0;
		System.out.println("Choose a method to run:");
		System.out.println("1. Exhaustive Search");
		System.out.println("2. Random Search");
		System.out.println("3. Genetic Algorithm Search");
		System.out.println("4. Genetic Algorithm Search");

		Scanner sc = new Scanner(System.in);
		choice = sc.nextInt();

		if (fileName.equals("data2.txt") && choice == 1) {

			int[] indexes = new int[labelsArray.length];

			for (int i = 0; i < labelsArray.length; i++) {

				indexes[i] = i;

			}// end i for

			long start = System.currentTimeMillis();
			printCombinations(indexes);
			long end = System.currentTimeMillis();

			exhaustTime = end - start;

			compute();

			System.out.println("The exhaust search time: " + exhaustTime / 1000
					+ " seconds");
			System.out.println("The total distances: " + totalDistances);
			System.out.println("The mean distances: " + meanDistances);
			System.out.println("The standard devition: " + stddev);

			System.out.println("Bins:");

			for (int i = 0; i < bins.length; i++) {

				System.out.println(i + " " + bins[i]);

			}// end i for

		}// end if for exhaust search

		if (fileName.equals("data2.txt") && choice == 2) {

			int n = 1000000;

			long start = System.currentTimeMillis();

			for (int i = 0; i < n; i++) {

				printRandomCombinations();

			}// end i for

			long end = System.currentTimeMillis();

			compute();

			randomSearchTime = end - start;

			System.out.println("The random search time: " + randomSearchTime
					+ " seconds");
			System.out.println("The total distances: " + totalDistances);
			System.out.println("The mean distances: " + meanDistances);
			System.out.println("The standard devition: " + stddev);

			System.out.println("Bins:");

			for (int i = 0; i < bins.length; i++) {

				System.out.println(i + " " + bins[i]);

			}// end i for

		}// end if for random search
		if (fileName.equals("data2.txt") && choice == 3) {

			geneticAlgorithm();

		}// end if for genetic algorithm

		if (fileName.equals("data2.txt") && choice == 4) {

			Permutation simTest = new Permutation();

			for (int i = 0; i < arrayPoints.length; i++) {

				simTest.setGene(i, arrayPoints[i].getLabel());

			}// end i for

			simulatedAnnealing(simTest);

		}// end if for genetic algorithm

	}// end run

	/**************************************************************************
	 * Creating the matrix
	 **************************************************************************/

	/**
	 * create the matrices
	 * 
	 * @return double[][] the specific matrix
	 */
	public double[][] createMatrix() {

		double[][] vector;
		Point2D point = new Point2D.Double(2.2, 2.2);
		int i = 0, j = 0, index = 0;
		char label = 0;
		double number = 0, number1 = 0, x = 0, y = 0;
		int row = 0, col = 0;
		double[][] createdMatrix = null;
		double temp[][] = null;
		String[] tempLabelsArray = null;

		System.out.println("Creating Matrix...");
		System.out
				.println("Enter the name of the matrix file: (include the type of file)");
		// get the file from the reader
		Scanner sc = new Scanner(System.in);
		fileName = sc.next();

		try {
			File file = new File(fileName);
			FileReader fin = new FileReader(file);
			BufferedReader matrixFile3 = new BufferedReader(fin);

			// Read the edges and insert
			String line;
			while ((line = matrixFile3.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(line);

				// count the number or columns
				col = st.countTokens();
				// counting the number of rows
				row++;

			}// end while

			// create a matrix with the counted row and col
			createdMatrix = new double[row][col];
			double[][] tempMatrix = new double[row - 1][col];

			System.out.println("Created a matrix of size : "
					+ createdMatrix.length + " by " + createdMatrix[0].length);

			// re-copy the matrix without the titles
			for (int k = 1; k < createdMatrix.length; k++) {

				for (int l = 0; l < createdMatrix[0].length; l++) {

					tempMatrix[k - 1][l] = createdMatrix[k][l];

				}// end l for

			}// end k for

			arrayPoints = new Points2DLabel[tempMatrix.length];
			labelsArray = new char[tempMatrix.length];

			System.out.println("Created a matrix of size : "
					+ tempMatrix.length + " by " + tempMatrix[0].length);

			fin = new FileReader(file);
			matrixFile3 = new BufferedReader(fin);

			// add in the numbers from the text file to the matrix
			while ((line = matrixFile3.readLine()) != null) {

				StringTokenizer st = new StringTokenizer(line);

				if (i != 0 && !line.equals("")) {

					if (fileName.equals("data.txt")
							|| fileName.equals("data2.txt")) {

						x = Double.parseDouble(st.nextToken());
						y = Double.parseDouble(st.nextToken());

						label = st.nextToken().charAt(0);
						labelsArray[i - 1] = label;

					}// end if

					point = new Point2D.Double(x, y);
					Points2DLabel pointLabel = new Points2DLabel(point, label,
							index);
					arrayPoints[i - 1] = pointLabel;

				}// end if

				index++;
				i++;

			}// end while

			System.out.println();

		} catch (IOException e) {

			System.err.println(e);

		}// end catch

		return temp;

	}// end createMatrix

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

	public void printCombinations(int[] str) throws IOException {

		V = str.length - 1;
		min = 3.29;
		max = 10.29;
		binWidth = (max - min) / numBins;

		// Permutation Tester

		val = new int[V + 1];
		for (int k = 0; k <= V; k++)
			val[k] = 0;

		p(0);

	}// end print combinations

	public void p(int k) {

		now++;
		val[k] = now;
		if (now == V)
			handleP();
		for (int i = 1; i <= V; i++)
			if (val[i] == 0)
				p(i);
		now--;
		val[k] = 0;

	}

	public void handleP() {

		count++;
		System.out.println("count" + count);

		printCost(val);

		for (int i = 0; i < val.length; i++) {

			System.out.print(val[i]);

		}

	}

	public void placeInBin(double perm) {

		int place = (int) ((perm - min) / binWidth);

		System.out.println(perm + " - " + min + " = " + (perm - min));

		System.out.println(place + " place");
		bins[place]++;
	}

	/**************************************************************************
	 * Random Search Algorithm Running - 1,000,000 random permutations
	 **************************************************************************/

	/**
	 * Finds 1,000,000 random permutations
	 */
	public void printRandomCombinations() {

		min = 3.29;
		max = 10.29;
		binWidth = (max - min) / numBins;

		Integer[] numbers = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };

		Collections.shuffle(Arrays.asList(numbers));

		int[] numbersTest = new int[numbers.length];

		for (int i = 0; i < numbers.length; i++) {

			numbersTest[i] = numbers[i];

		}// end i for

		for (int i = 0; i < numbers.length; i++) {

			System.out.print(numbersTest[i] + " ");

		}

		printCost(numbersTest);

	}// end printRandomCombinations

	/**************************************************************************
	 * Write the costs of each permutation to a text file
	 **************************************************************************/

	/**
	 * Takes the file that contain the permutation Computers each cost of travel
	 * for each permutation
	 * 
	 * @throws IOException
	 */
	public void printCost(int[] s) {

		double totalDistanceTraveled = 0.0;
		double distance = 0.0;

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

		totalDistances += totalDistanceTraveled;
		countDistances++;

		if (totalDistanceTraveled > max)

			max = totalDistanceTraveled;

		if (totalDistanceTraveled < min)

			min = totalDistanceTraveled;

		xisq = Math.pow(totalDistanceTraveled, 2);
		xisqSum += xisq;

		xiSumsq += totalDistanceTraveled;

		placeInBin(totalDistanceTraveled);
		// out = new PrintWriter(new FileWriter("distanceCost.txt", true));
		// out.println(totalDistanceTraveled + ",");
		// out.close();
		// out.flush();
		totalDistanceTraveled = 0.0;

	}// end printCost

	/**************************************************************************
	 * Genetic Algorithm - Methods: computeFitnessGenetic(Permutation [] a) ,
	 * findBestIndividual(Permutation[] a) , tournamentSlect() ,
	 * crossover(Permutation indiv1, Permutation indiv2) , evaluate(Permutation
	 * [] a) , setNewPopulation(Permutation [] a). Steps: 1. Create a random
	 * list of individuals for a population of any given size 2. Compute the top
	 * 5 individuals of the previous population and save it to the new
	 * population 3. Fill in the new population with individuals a. Pick the
	 * best individual for individual 1 amongst a tournament selection of
	 * individuals from the original population b. Pick the best individual for
	 * individual 1 amongst a tournament selection of individuals from the
	 * original population c. If the random generator number falls < less the
	 * rates (crossover, both mutations) perform the given method d. Add
	 * individual 1 & 2 to the new population 4. Increment the counter to point
	 * to the next empty space in the new population 5. Repeat the steps until
	 * all of the new population spaces have been filled 6. Calculate the
	 * fittest individual from the new population 7. Add to a array of best
	 * individuals 8. Repeat Steps 1 - 7 for the number of tests while saving
	 * the best 5 individuals from the previous population 9. After tests are
	 * complete, find the best individual from the array of best individuals 10.
	 * Genetic Algorithm Complete
	 **************************************************************************/

	/**************************************************************************
	 * Compute the distances of each permutation
	 **************************************************************************/

	public void compute() {

		stddev = 0.0;

		meanDistances = totalDistances / countDistances;

		// compute standard dev.

		xiSumsq = Math.pow(xiSumsq, 2);

		stddev = Math.sqrt((xisqSum - (xiSumsq) / countDistances)
				/ (countDistances - 1));

	}// end getValues

	/**
	 * Method to run the genetic algorithm
	 */
	public void geneticAlgorithm() {

		min = 3.29;
		max = 10.29;
		binWidth = (max - min) / numBins;

		topIndividuals = 5;
		populationSize = 50 + topIndividuals; // population size
		numberOfTests = 2000; // max number of iterations
		mutationRate = 0.05; // probability of mutation
		crossoverRate = 0.7; // probability of crossover

		int count; // to keep count of number of individuals in a given
					// population
		int countBestIndiv = 0;
		Permutation[] newPop = new Permutation[populationSize];
		Permutation[] indiv = new Permutation[2];
		bestIndividuals = new Permutation[2000];// find the best 2000
												// individuals
		// create the population
		mainPopulation = new Permutation[populationSize];

		// initial population
		// fill the population with random individuals
		for (int i = 0; i < populationSize; i++) {

			mainPopulation[i] = new Permutation();
			mainPopulation[i].randGenes(labelsArray);

		}// end i for

		double totalFit = evaluate(mainPopulation);

		// current population
		System.out.print("Total Fitness = " + totalFit);
		System.out.println(" ; Best Fitness = "
				+ findBestIndividual(mainPopulation).getFitnessValue());

		for (int iter = 0; iter < numberOfTests; iter++) {
			count = 0;

			binWidth = (max - min) / numBins;

			// /Start the genetic algorithm

			printPop(mainPopulation);
			// find the top 10% individuals of the new population and keep them
			computeFitnessGenetic(mainPopulation);

			// fill the first 5 spots of the new population with the best 5
			// individuals from the last population
			for (int i = topIndividuals; i > 0; i--) {

				newPop[count] = top5[count];

				count++;
			}

			// init population
			for (int i = count; i < populationSize; i++) {

				newPop[i] = new Permutation();
				mainPopulation[i].randGenes(labelsArray);

			}// end i for

			while (count < populationSize) {

				// Selection method
				indiv[0] = tournamentSelection();
				indiv[1] = tournamentSelection();

				// Crossover method
				if (m_rand.nextDouble() < crossoverRate) {

					indiv = crossover(indiv[0], indiv[1]);

				}

				// Mutation method
				if (m_rand.nextDouble() < mutationRate) {

					indiv[0].mutate();
				}

				if (m_rand.nextDouble() < mutationRate) {

					indiv[1].mutate();
				}

				// add to the new population
				newPop[count] = indiv[0];
				newPop[count + 1] = indiv[1];
				count += 2;

			}// end while

			// set the new population to the current population
			setNewPopulation(newPop);

			// reevaluate current population
			totalfitness = evaluate(mainPopulation);
			Permutation bestindiv = findBestIndividual(newPop);

			System.out.print("Total Fitness = " + totalfitness);
			System.out.println(" ; Best Fitness = "
					+ bestindiv.getFitnessValue());

			// add the best individual in current population to array of best
			// individuals
			bestIndividuals[countBestIndiv] = bestindiv;
			countBestIndiv++;

			System.out.println("best individual" + bestindiv.getFitnessValue());
			System.out.println("Max: " + max);
			System.out.println("Min: " + min);

			if (bestindiv.getFitnessValue() > max) {

				max = bestindiv.getFitnessValue();

				System.out.print(max + " changed");
			}

			if (bestindiv.getFitnessValue() < min) {

				min = bestindiv.getFitnessValue();

				System.out.print(min + " changed");
			}

			placeInBin(bestindiv.getFitnessValue());

		}// end for

		Permutation leader = findBestIndividual(bestIndividuals);

		System.out.println("\nThe leader of the pack is: "
				+ leader.chartoString() + " with a fitness value of: "
				+ leader.getFitnessValue());

	}// end genetic algorithm

	/**
	 * Method to compute the top 5 best individuals and save them for the next
	 * population
	 * 
	 * @param a
	 *            the population
	 */
	public void computeFitnessGenetic(Permutation[] a) {

		System.out.println("Compute Fitness Genetic");

		fitnesses = new double[a.length];
		top5 = new Permutation[5];
		double[] top5fitness = new double[5];

		evaluate(a);

		// get the fitness values for the current population
		for (int i = 0; i < a.length; i++) {

			System.out.println(a.length);
			System.out.println(fitnesses[i] + "-->" + a[i].getFitnessValue());
			fitnesses[i] = a[i].getFitnessValue();

		}// end i for

		// sort the fitness values
		Arrays.sort(fitnesses);

		int j = a.length - 5;
		int m = 0;

		// get the top 5 fitness value from the array
		while (j < a.length) {

			top5fitness[m] = fitnesses[j];
			j++;
			m++;

		}// end while

		// search through the population to get the individual
		// set the index in top5 array where the individual have been found in
		// the original population
		for (int n = 0; n < a.length; n++) {

			if (top5fitness[0] == a[n].getFitnessValue()) {

				top5[0] = a[n];

			}// end if
			if (top5fitness[1] == a[n].getFitnessValue()) {

				top5[1] = a[n];

			}// end if
			if (top5fitness[2] == a[n].getFitnessValue()) {

				top5[2] = a[n];

			}// end if
			if (top5fitness[3] == a[n].getFitnessValue()) {

				top5[3] = a[n];

			}// end if
			if (top5fitness[4] == a[n].getFitnessValue()) {

				top5[4] = a[n];

			}// end if

		}// end n for

		System.out.println("Print top 5");
		printPop(top5);
		System.out.println("Print top 5");

	}// end computer fitness

	/**
	 * Get the best individual
	 * 
	 * @return
	 */
	public Permutation findBestIndividual(Permutation[] a) {

		int max = 0, min = 0;
		double currentMax = 0.0;
		double currentMin = 1.0;
		double currentValue;

		for (int i = 0; i < a.length; i++) {
			currentValue = a[i].getFitnessValue();
			if (currentMax < currentMin) {

				currentMax = currentMin = currentValue;
				max = min = i;

			}// end if
			if (currentValue > currentMax) {

				currentMax = currentValue;
				max = i;

			}// end if

			if (currentValue < currentMin) {

				currentMin = currentValue;
				min = i;

			}// end if
		}// end idx for

		return a[min]; // minimization

	}

	/**
	 * Selects the best individuals for the population
	 * 
	 * @return a individual
	 */
	public Permutation tournamentSelection() {
		// Create a tournament population
		Permutation[] tournament = new Permutation[5];
		// For each place in the tournament get a random individual
		for (int i = 0; i < tournament.length; i++) {

			int randomId = (int) (Math.random() * populationSize);

			tournament[i] = mainPopulation[randomId];

		}

		// evaluate the tournament population
		double totalFit = evaluate(tournament);

		// current population
		System.out.print("Total Fitness = " + totalFit);
		System.out.println(" ; Best Fitness = "
				+ findBestIndividual(tournament).getFitnessValue());

		// Get the fittest
		Permutation fittest = findBestIndividual(tournament);

		System.out.println("Fittest in tournament: " + fittest.chartoString());
		System.out.println("Fittest Value: " + fittest.getFitnessValue());

		return fittest;

	}

	/**
	 * Performs uniform crossover between two individuals
	 * 
	 * @param indiv1
	 *            the first selected indidual
	 * @param indiv2
	 *            the second selected individual
	 * @return
	 */
	public Permutation[] crossover(Permutation individual1,
			Permutation individual2) {

		Permutation[] newIndividual = new Permutation[2];

		newIndividual[0] = new Permutation();
		newIndividual[1] = new Permutation();

		// loop through the genes
		for (int i = 0; i < individual1.length(); i++) {

			// performing crossover
			if (Math.random() <= uniformRate) {

				newIndividual[0].setGene(i, individual1.getGene(i));
				newIndividual[1].setGene(i, individual2.getGene(i));

			} else {

				newIndividual[0].setGene(i, individual2.getGene(i));
				newIndividual[1].setGene(i, individual1.getGene(i));
			}
		}// end i for

		return newIndividual;

	}// end crossover

	/**
	 * Evaluate the current population
	 * 
	 * @return
	 */
	public double evaluate(Permutation[] a) {

		double totalfitness1 = 0.0;

		for (int i = 0; i < a.length; i++) {

			totalfitness1 += a[i].evaluate(arrayPoints);

		}// end i for

		return totalfitness1;

	}// end evaluate

	/**
	 * Set the new population to the default population Used to evaluate the
	 * fitness of the new population
	 * 
	 * @param newPop
	 */
	public void setNewPopulation(Permutation[] newPop) {

		// this.mainPopulation = newPop;

		for (int i = 0; i < newPop.length; i++) {

			mainPopulation[i] = newPop[i];

		}// end i for

	}// end setPopulation

	/**************************************************************************
	 * Simulating Annealing
	 **************************************************************************/

	/**
	 * Compute the simulated annealing
	 * 
	 * @param a
	 *            the current permutation
	 * @throws IOException
	 */
	public void simulatedAnnealing(Permutation a) throws IOException {

		min = 3.29;
		max = 10.29;
		binWidth = (max - min) / numBins;

		int iteration = -1;

		double temperature = 100.0;
		double deltaDistance = 0;
		double coolingRate = 0.99;
		double absoluteTemperature = 0.00001;
		double random = 0;

		double distance = a.evaluate(arrayPoints);
		// System.out.println("The distance is: "+distance);

		while (temperature > absoluteTemperature) {

			Permutation nextOrder = randomSwitch(a);

			deltaDistance = nextOrder.evaluate(arrayPoints) - distance;

			// if the new order has a smaller distance
			// or if the new order has a larger distance but
			// satisfies Boltzman condition then accept the arrangement

			random = Math.random();

			if ((deltaDistance < 0)
					|| (distance > 0 && Math.exp(-deltaDistance / temperature) > random)) {
				for (int i = 0; i < nextOrder.length(); i++)

					a.setGene(i, nextOrder.getGene(i));

				distance = deltaDistance + distance;

				if (distance > max)
					max = distance;

				if (distance < min)
					min = distance;

				placeInBin(distance);

			}

			// cool down the temperature
			temperature *= coolingRate;

			iteration++;
		}

		System.out.println("The shortest distance is " + distance);
		System.out.println("Number of iterations " + iteration);

	}// end simulated annealing

	/**
	 * Does random switching
	 * 
	 * @param a
	 * @return
	 */
	public Permutation randomSwitch(Permutation a) {

		Permutation SWIT = new Permutation();

		int rand1 = (int) (Math.random() * a.length());
		int rand2 = (int) (Math.random() * a.length());

		char temp = a.getGene(rand1);
		a.setGene(rand1, a.getGene(rand2));
		a.setGene(rand2, temp);

		for (int z = 0; z < a.length(); z++) {

			SWIT.setGene(z, a.getGene(z));

		}// end z for

		return SWIT;

	}// end random switch

	/**************************************************************************
	 * Computations and Manipulations
	 **************************************************************************/

	/**
	 * Chages the String[] to a String
	 * 
	 * @param a
	 *            the String[]
	 * @return
	 */
	public String changeToString(char[] a) {

		String concat = "";

		StringBuffer sb = new StringBuffer(concat);

		for (int i = 0; i < a.length; i++) {

			sb.append(a[i]);

		}// end i for

		concat = sb.toString();

		return concat;

	}// end changeToChar

	/**
	 * compute the mean of the distances
	 * 
	 * @param a
	 *            the array of distances
	 * @return
	 */
	public double computeMeanDistances(double[] a) {

		double totalDistances = 0.0;
		double meanDistances = 0.0;

		for (int i = 0; i < a.length; i++) {

			totalDistances += a[i];

		}// end i for

		meanDistances = totalDistances / a.length;

		return meanDistances;

	}// end compute mean distances

	/**
	 * computes the standard deviation of the array of distances
	 * 
	 * @return
	 */
	public double computeStDDev(double[] a) {

		double xi = 0.0;
		double xisq = 0.0;
		double xisqSum = 0.0;
		double xiSumsq = 0.0;
		double stddev = 0.0;
		int n = a.length;

		for (int i = 0; i < a.length; i++) {

			xisq = Math.pow(a[i], 2);
			xisqSum += xisq;
			xiSumsq += a[i];

		}// end i for

		xiSumsq = Math.pow(xiSumsq, 2);
		stddev = Math.sqrt((xisqSum - (xiSumsq) / n) / (n - 1));

		return stddev;

	}// end computerStDDev

	/**
	 * Print the population
	 * 
	 * @param a
	 *            the population
	 */
	public void printPop(Permutation[] a) {

		for (int i = 0; i < a.length; i++) {

			System.out.println(a[i].chartoString());

		}// end i for
	}// end printPop
}// end class
