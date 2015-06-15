import java.awt.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.StringTokenizer;

public class Permutation {

	private String genesString;
	private double fitnessValue;
	static int SIZE = 6;
	private static final double mutationRate = 0.015;
	private char[] genes = new char[14];
	private double test = 1;

	public Permutation() {

	}// end constructor

	public double getFitnessValue() {

		return this.fitnessValue;
	}

	public void setFitnessValue(int fitnessValue) {

		this.fitnessValue = fitnessValue;
	}

	public char getGene(int index) {

		return this.genes[index];

	}

	public char[] returnGene() {

		return this.genes;

	}

	public int length() {

		return genes.length;

	}// end size

	public void setGene(int index, char gene) {

		this.genes[index] = gene;

		// fitnessValue = 0.0;

	}

	public void randGenes(char[] a) {

		this.genes = new char[a.length];

		java.util.List<char[]> asList = Arrays.asList(a);
		ArrayList<Character> listC = new ArrayList<Character>();

		for (char c : a) {

			listC.add(c);

		}

		Collections.shuffle(listC);

		Character[] array = listC.toArray(new Character[listC.size()]);

		for (int i = 0; i < array.length; i++) {

			this.setGene(i, array[i]);

		}// end i for
	}// end ran genes

	public void mutate() {

		// Loop through genes

		Random rand = new Random();

		int index = rand.nextInt(genes.length);
		int index2 = rand.nextInt(genes.length);

		while (index == index2) {

			index = rand.nextInt(genes.length);
		}

		char temp = this.getGene(index);
		this.setGene(index, this.genes[index2]);
		this.setGene(index2, temp);

	}// end mutate

	public String chartoString() {

		this.genesString = "";

		StringBuffer sb = new StringBuffer(genesString);

		for (int i = 0; i < this.genes.length; i++) {

			sb.append(this.genes[i]);

		}// end i for

		this.genesString = sb.toString();

		return genesString;

	}// end charToString

	public double evaluate(Points2DLabel[] arrayPoints) {

		double fitness = 0.0;
		double distance = 0.0;

		this.genesString = chartoString();

		String currentPermutation = genesString;

		StringTokenizer st = new StringTokenizer(currentPermutation);

		StringBuffer sb = new StringBuffer(currentPermutation);
		Points2DLabel currentPoint2D = null;
		Points2DLabel previousPoint2D = null;

		for (int j = 0; j < sb.length(); j++) {

			char currentLabelList = sb.charAt(j);

			for (int k = 0; k < arrayPoints.length; k++) {

				if (currentLabelList == arrayPoints[k].getLabel()) {

					if (previousPoint2D == null) {

						previousPoint2D = arrayPoints[k];

					} else {

						currentPoint2D = arrayPoints[k];

						// calculate the distance between the previous
						// points and the current point
						distance = previousPoint2D.getPoint().distance(
								currentPoint2D.getPoint());

						fitness += distance;

						previousPoint2D = currentPoint2D;

					}// end else

				}// end if

			}// end k for

			distance = previousPoint2D.getPoint().distance(
					arrayPoints[0].getPoint());

			fitness += distance;

		}// end j for

		this.fitnessValue = fitness;
		return fitnessValue;
	}

}// end class
