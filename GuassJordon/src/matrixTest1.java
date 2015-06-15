import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.StringTokenizer;

/**
 * Test 1 Spring 2014 Taglirini When the program asks you to enter a text file.
 * 1. Type "data.txt" Then the program will ask you to type in another file This
 * file is the systems of equation for the second part 2. Type "system.txt"
 * 
 * @author Kavy Rattana
 * 
 */
public class matrixTest1 {

	public int choice, constantLength, inverConLength, countTransposeRow,
			countTransposeCol, classifiedCountC1, classifiedCountC2,
			misclassifiedCountC1, misclassifiedCountC2;
	public String classified, classifiedAs;
	public double scalar, determinantNum, r, discriminantG1, discriminantG2,
			determinantM1, determinantM2, episilon;
	public double[][] inverseM1, inverseM2, idenMatrixM1, idenMatrixM2,
			augmentedMatrix, productMatrix, sumMatrixM1, sumMatrixM2,
			determinantMatrixM1, determinantMatrixM2, temp, mean1, mean2;
	public double[][] vector, vector2, tempVectorV1, tempVectorV2,
			transposeMatrix, covarianceMatrixM1, covarianceMatrixM2,
			identityMatrix, augMatrixM1, augMatrixM2, idenTemp;
	public boolean matrix1Created, class1, class2;
	public Scanner input2;
	public ArrayList<double[][]> arrays = new ArrayList<double[][]>(),
			arrays2 = new ArrayList<double[][]>(),
			arrays3 = new ArrayList<double[][]>();
	public ArrayList<Double> discriminantArrayListM1 = new ArrayList<Double>(),
			discriminantArrayListM2 = new ArrayList<Double>();
	public double[][] test;
	public ArrayList<Integer> testInts;

	/**
	 * Main to initialize the class
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

		matrixTest1 obj = new matrixTest1();
		obj.run();

	}// end main

	/**
	 * Method to run the program
	 */
	public void run() {

		// initialize the variables
		double[][] temp = createMatrix();
		sumMatrixM1 = new double[2][2];
		sumMatrixM2 = new double[2][2];
		covarianceMatrixM1 = new double[2][2];
		covarianceMatrixM2 = new double[2][2];
		determinantM1 = 0.0;
		determinantM2 = 0.0;
		double[][] tempCovM1 = new double[2][2];
		double[][] tempCovM2 = new double[2][2];
		double[][] tempTrans = new double[1][2];

		computeMean(temp);

		double[][] subTemp1 = new double[2][1];

		// calculate the covariance of M1
		for (int i = 0; i < arrays.size(); i++) {

			// get the matrix from the arraylist of matrices
			tempVectorV1 = arrays.get(i);
			// subtract the mean vector from every vector in class 2
			subTemp1 = subtraction(tempVectorV1, transpose(mean1));
			// compute the product of the vector and it's transpose
			tempTrans = multiply(subTemp1, transpose(subTemp1));
			// add the matrix to sum of matrices
			addMatrixM1(tempTrans);

		}// end i for

		double[][] subTemp2 = new double[2][1];

		// calculate the covariance of M2
		for (int i = 0; i < arrays2.size(); i++) {
			// get the matrix from the arraylist of matrices

			tempVectorV2 = arrays2.get(i);
			// subtract the mean vector from every vector in class 2
			subTemp2 = subtraction(tempVectorV2, transpose(mean2));
			// compute the product of the vector and it's transpose
			tempTrans = multiply(subTemp2, transpose(subTemp2));

			// add the matrix to sum of matrices
			addMatrixM2(tempTrans);

		}// end i for

		// multiply the sum matrix M1 by the scalar
		System.out.println("The sum of M1");
		debug(sumMatrixM1);
		System.out.println("The sum of M2");
		debug(sumMatrixM2);
		System.out.println();
		System.out.println("The scalar is: " + scalar);
		covarianceMatrixM1 = scalarMultiplication(sumMatrixM1, scalar);
		// multiply the sum matrix M2 by the scalar
		covarianceMatrixM2 = scalarMultiplication(sumMatrixM2, scalar);

		// Print out the covariance matrices
		System.out.println("\n\nClass 1\t\tSigma 1");
		debug(covarianceMatrixM1);
		System.out.println("Class 2\t\tSigma 2");
		debug(covarianceMatrixM2);

		// copy the covariance to be used again
		tempCovM1 = copyMatrix(covarianceMatrixM1);
		tempCovM2 = copyMatrix(covarianceMatrixM2);

		System.out.println();
		System.out.println();

		// solve for the determinant using gaussian implementation
		determinantMatrixM1 = determinant(covarianceMatrixM1);
		determinantMatrixM2 = determinant(covarianceMatrixM2);

		// print the determinants of the classes
		determinantM1 = findDelta(determinantMatrixM1, r);
		determinantM2 = findDelta(determinantMatrixM2, r);
		System.out.println("\n\n|Sigma1|= " + determinantM1 + "\t\t|Sigma2|= "
				+ determinantM2 + "\n\n");

		// make an identity matrix for the inputed matrix
		idenMatrixM1 = makeIden(tempCovM1);
		idenMatrixM2 = makeIden(tempCovM2);

		// combines the inputed matrix with the identity matrix
		augMatrixM1 = createAugmentedMatrix(tempCovM1, idenMatrixM1);
		augMatrixM2 = createAugmentedMatrix(tempCovM2, idenMatrixM2);

		// compute the inverse of the augmented matrices
		inverseM1 = partitionInverse(inverse(augMatrixM1));
		inverseM2 = partitionInverse(inverse(augMatrixM2));

		// prints the output of the inverse
		System.out.println("\n\nSigma1Inv");
		printMatrix(inverseM1);
		System.out.println("\nSigma2Inv");
		printMatrix(inverseM2);

		System.out.println("g for m1 and m2");
		System.out.println("\tx\t\ty\tClass\tClassified as\t\tg1\t\t\tg2");

		// for M1
		debug(mean1);
		discriminantG1 = discriminant1(mean1, inverseM1, transpose(mean1),
				determinantM1);
		discriminantG2 = discriminant1(mean2, inverseM2, transpose(mean1),
				determinantM2);

		if (discriminantG1 > discriminantG2) {

			classifiedAs = "C1";
			System.out.println(mean1[0][0] + "\t" + mean1[0][1] + "\tC1\t"
					+ classifiedAs + "\t" + discriminantG1 + "\t"
					+ discriminantG2);

		} else if (discriminantG1 < discriminantG2) {

			classifiedAs = "C2";
			System.out.println(mean1[0][0] + "\t" + mean1[0][1] + "\tC1\t"
					+ classifiedAs + "\t" + discriminantG1 + "\t"
					+ discriminantG2);

		}// end else if

		// for M2
		debug(mean2);
		discriminantG1 = discriminant1(mean1, inverseM1, transpose(mean2),
				determinantM1);
		discriminantG2 = discriminant1(mean2, inverseM2, transpose(mean2),
				determinantM2);

		if (discriminantG1 > discriminantG2) {

			classifiedAs = "C1";
			System.out.println(mean2[0][0] + "\t" + mean2[0][1] + "\tC1\t"
					+ classifiedAs + "\t" + discriminantG1 + "\t"
					+ discriminantG2);

		} else if (discriminantG1 < discriminantG2) {

			classifiedAs = "C2";
			System.out.println(mean2[0][0] + "\t" + mean2[0][1] + "\tC1\t"
					+ classifiedAs + "\t" + discriminantG1 + "\t"
					+ discriminantG2);

		}// end else if

		System.out.println("\n\nPoint");

		System.out.println("Class 1");
		System.out.println("\tx\t\ty\tClass\tClassified as\t\tg1\t\t\tg2");
		for (int i = 0; i < arrays.size(); i++) {

			tempVectorV1 = arrays.get(i);

			class1 = true;
			// g1 function
			discriminantG1 = discriminant1(mean1, inverseM1, tempVectorV1,
					determinantM1);
			// g2 function
			discriminantG2 = discriminant1(mean2, inverseM2, tempVectorV1,
					determinantM2);

			// checks to see if g1 > g2 then it is in the current class
			// if not, then it is in the other class
			if (discriminantG1 > discriminantG2 && class1 == true) {

				classified = "Classified";
				classifiedAs = "C1";
				System.out.println(tempVectorV1[0][0] + "\t"
						+ tempVectorV1[1][0] + "\tC1\t" + classifiedAs + "\t"
						+ discriminantG1 + "\t" + discriminantG2 + "\t"
						+ classified);
				classifiedCountC1++;

				// checks to see if g1 > g2 then it is in the current class
				// if not, then it is in the other class
			} else if (discriminantG1 < discriminantG2 && class1 == true) {

				classified = "Misclassified";
				classifiedAs = "C2";
				System.out.println(tempVectorV1[0][0] + "\t"
						+ tempVectorV1[1][0] + "\tC1\t" + classifiedAs + "\t"
						+ discriminantG1 + "\t" + discriminantG2 + "\t"
						+ classified);
				misclassifiedCountC1++;

			}// end else if

			class1 = false;

		}// end i for

		System.out.println("\nClass 2");
		System.out.println("\tx\t\ty\tClass\tClassified as\t\tg1\t\t\tg2");

		for (int j = 0; j < arrays2.size(); j++) {

			// get the vector from the arraylist
			tempVectorV2 = arrays2.get(j);

			class2 = true;

			// g1 function
			discriminantG1 = discriminant1(mean1, inverseM1, tempVectorV2,
					determinantM1);
			// g2 function
			discriminantG2 = discriminant1(mean2, inverseM2, tempVectorV2,
					determinantM2);

			// checks to see if g1 > g2 then it is in the current class
			// if not, then it is in the other class
			if (discriminantG1 < discriminantG2 && class2 == true) {

				classified = "Classified";
				classifiedAs = "C2";
				System.out.println(tempVectorV2[0][0] + "\t"
						+ tempVectorV2[1][0] + "\tC2\t" + classifiedAs + "\t"
						+ discriminantG1 + "\t" + discriminantG2 + "\t"
						+ classified);
				classifiedCountC2++;

				// checks to see if g1 > g2 then it is in the current class
				// if not, then it is in the other class
			} else if (discriminantG1 > discriminantG2 && class2 == true) {

				classified = "Misclassified";
				classifiedAs = "C1";
				System.out.println(tempVectorV2[0][0] + "\t"
						+ tempVectorV2[1][0] + "\tC2\t" + classifiedAs + "\t"
						+ discriminantG1 + "\t" + discriminantG2 + "\t"
						+ classified);
				misclassifiedCountC2++;

			}// end else if

			class2 = false;

		}// end j for

		// print out the classified and misclassified points for class 1 and
		// class 2
		System.out.println("The # of correct classified in C1: "
				+ classifiedCountC1);
		System.out.println("The # of misclassified in C1: "
				+ misclassifiedCountC1);
		System.out.println("The # of correct classified in C2: "
				+ classifiedCountC2);
		System.out.println("The # of  misclassified in C2: "
				+ misclassifiedCountC2);

		// contour boundary
		contour(mean1, inverseM1, determinantM1, mean2, inverseM2,
				determinantM2);

		double[][] vecTemp = new double[2][1];
		System.out.println(arrays3.size());

		// print out the contour boundary points
		for (int i = 0; i < arrays3.size(); i++) {

			vecTemp = arrays3.get(i);
			System.out.print(vecTemp[0][0] + "\t" + vecTemp[1][0]);
			System.out.println();

		}// end i for

		// create copies of the Coefficient A matrix
		double[][] systemsMatrix = createMatrix();
		double[][] copyMatrix1 = partitionSystems(systemsMatrix);
		double[][] copyMatrix2 = partitionSystems(systemsMatrix);
		double[][] copyMatrix3 = partitionSystems(systemsMatrix);
		double[][] copyMatrix4 = partitionSystems(systemsMatrix);
		double[][] copyMatrix5 = partitionSystems(systemsMatrix);
		double[][] copyMatrix6 = partitionSystems(systemsMatrix);

		// guassJordan of system A
		System.out.println("\nGuass Jordan Row Reduction of A");
		printOutput(guassJordan(systemsMatrix));
		// determinant of coefficient A
		double deterSystems = findDelta(determinant(copyMatrix1), r);
		System.out.println("\nDeterminant of Coefficients A: " + deterSystems);

		// inverse of coefficient A
		System.out.println("\nInverse of Coefficient A");
		double[][] aInverse = partitionInverse(inverse(createAugmentedMatrix(
				copyMatrix2, makeIden(copyMatrix2))));
		debug(aInverse);

		// determinant of A inverse
		System.out.println("\nDeterminant of Coefficient A Inverse");
		double deterSystemsInver = findDelta(
				determinant(partitionInverse(inverse(createAugmentedMatrix(
						copyMatrix3, makeIden(copyMatrix3))))), r);

		System.out.println("\nDeterminant of Systems Inverse: "
				+ deterSystemsInver);

		// multiply the determinant of Coefficient A with its inverse
		System.out.println("deterSystems: " + deterSystems
				+ "\ndeterSystemsInver: " + deterSystemsInver);
		double productDeter = deterSystems * deterSystemsInver;
		System.out
				.println("\nDeterminant A * Determinant Coefficient A Inverse: "
						+ productDeter);

		// checking if the Coefficient A Inverse exist
		System.out.println("\nChecking Coefficient A Inverse\n");
		System.out
				.println("\nChecking the Coefficient A Inverse with the Coefficients");
		printMatrix(multiply(copyMatrix4, aInverse));

		// computer the condition number
		System.out.println("\nComputing the Condition Number: ");
		System.out.println(conditionNumber(copyMatrix5, aInverse));

	}// end run

	/**
	 * create the matrices
	 * 
	 * @return double[][] the specific matrix
	 */
	public double[][] createMatrix() {

		System.out.println("Creating Matrix...");
		System.out
				.println("Enter the name of the matrix file: (include the type of file)");

		int i = 0, j = 0;

		// get the file from the reader
		Scanner sc = new Scanner(System.in);
		String fileName = sc.next();

		double number = 0;
		int row = 0, col = 0;

		// initialize the double [][] arrays
		double[][] createdMatrix = null;
		double temp[][] = null;

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

			System.out.println("Created a matrix of size : " + row + " by "
					+ col);

			fin = new FileReader(file);
			matrixFile3 = new BufferedReader(fin);

			// add in the numbers from the text file to the matrix
			while ((line = matrixFile3.readLine()) != null) {

				StringTokenizer st = new StringTokenizer(line);

				for (j = 0; j < col; j++) {

					if (i != 0) {

						if (fileName.equals("data.txt")
								|| fileName.equals("data2.txt")) {

							number = Double.parseDouble(st.nextToken());

						} else if (fileName.equals("system.txt")) {

							number = Integer.parseInt(st.nextToken());

						}// end else if

						createdMatrix[i][j] = number;

					}// end if

				}// end j for

				i++;

				System.out.println();

			}// end while

			temp = new double[row - 1][col];

			// re-copy the matrix without the titles
			int m = 1;
			for (int k = 1; k < createdMatrix.length; k++) {

				for (int l = 0; l < createdMatrix[0].length; l++) {

					temp[k - 1][l] = createdMatrix[k][l];

				}// end l for

			}// end k for

			if (fileName.equals("data.txt") || fileName.equals("data2.txt")) {

				vector = new double[temp[0].length / 2][1];
				vector2 = new double[temp[0].length / 2][1];

				double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;

				// make arrays for the vectors
				for (int n = 0; n < temp.length; n++) {

					vector = new double[temp[0].length / 2][1];
					vector2 = new double[temp[0].length / 2][1];

					for (int o = 0; o < temp[0].length; o++) {

						if (o == 0) {

							x1 = temp[n][o];

						}// end if
						if (o == 1) {

							y1 = temp[n][o];

						}// end if

						if (o == 2) {

							x2 = temp[n][o];

						}// end if
						if (o == 3) {

							y2 = temp[n][o];

						}// end if

					}// end o for

					vector[0][0] = x1;
					vector[1][0] = y1;
					vector2[0][0] = x2;
					vector2[1][0] = y2;

					arrays.add(vector);
					arrays2.add(vector2);

				}// end n for

				System.out.println("Class 1");
				for (double[][] c : arrays) {

					for (int p = 0; p < c.length; p++) {

						for (int q = 0; q < c[0].length; q++) {

							System.out.print("[" + c[p][q] + "]");

						}// end q
					}// end p

					System.out.println();

				}// end for
				System.out.println();

				System.out.println("Class 2");
				for (double[][] c : arrays2) {

					for (int p = 0; p < c.length; p++) {

						for (int q = 0; q < c[0].length; q++) {

							System.out.print("[" + c[p][q] + "]");

						}// end q
					}// end p

					System.out.println();
				}// end for

				System.out.println("Array 1 size: " + arrays.size());
				System.out.println("Array 2 size: " + arrays2.size());
			}// end if

		} catch (IOException e) {
			System.err.println(e);
		}// end catch

		return temp;

	}// end createMatrix

	/**
	 * prints out the matrices
	 * 
	 * @param a
	 *            the matrix to printr
	 */
	public void printMatrix(double[][] a) {

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				System.out.print(a[i][j] + "\t");

			}// end j for

			System.out.println();

		}// end i for

		System.out.println(a.length + " by " + a[0].length + " Matrix");

	}// end print

	/**
	 * Debugs the current matrix
	 * 
	 * @param d
	 *            the current matrix
	 */
	public void debug(double[][] d) {

		for (int i = 0; i < d.length; i++) {

			for (int k = 0; k < d[0].length; k++) {

				System.out.print(d[i][k] + "\t");
			}// end k for

			System.out.println();
		}// end i for

		System.out.println("Dimensions:" + "\t(" + d.length + " x "
				+ d[0].length + ")");

	}// end debug

	/**
	 * Computers the mean for given matrices.
	 * 
	 * @param a
	 *            the current matrix being addded
	 */
	public void computeMean(double[][] a) {

		double[][] A = a;
		double sumX1 = 0, sumY1 = 0, sumX2 = 0, sumY2 = 0, multiple;
		scalar = Math.pow(A.length - 1, -1);
		System.out.println(A.length - 1);
		System.out.println("Scalar: " + scalar);

		mean1 = new double[1][A[0].length / 2];
		mean2 = new double[1][A[0].length / 2];

		// mean for x1
		for (int i = 0; i < A.length; i++) {

			for (int j = 0; j < A[0].length; j++) {

				// if statements to sum up each column
				if (j == 0) {

					sumX1 += A[i][j];

				}// end j

				if (j == 1) {

					sumY1 += A[i][j];

				}// end j

				if (j == 2) {

					sumX2 += A[i][j];

				}// end j

				if (j == 3) {

					sumY2 += A[i][j];

				}// end j

			}// end i for

		}// end computer mean

		// multiply the mean by the scalar
		mean1[0][0] = sumX1 * scalar;
		mean1[0][1] = sumY1 * scalar;

		System.out.println("\nMU1");
		System.out.println("mx1\t\t\tmy1");

		// print out the mean1
		for (int k = 0; k < mean1.length; k++) {

			for (int l = 0; l < mean1[0].length; l++) {

				System.out.print(mean1[k][l] + "\t");

			}// end l for

		}// end k for

		if (A[0].length >= 4) {
			// multiply the mean by the scalar
			mean2[0][0] = sumX2 * scalar;
			mean2[0][1] = sumY2 * scalar;

			System.out.println("\n\nMU2");
			System.out.println("mx2\t\t\tmy2");

			// print out mean2
			for (int k = 0; k < mean2.length; k++) {

				for (int l = 0; l < mean2[0].length; l++) {

					System.out.print(mean2[k][l] + "\t");

				}// end l for

			}// end k for
		}// end if
	}// end compute mean

	/**
	 * Computers the transpose for a matrix
	 * 
	 * @param a
	 *            the matrix
	 * @return the transposed matrix
	 */
	public double[][] transpose(double[][] a) {

		double[][] A = a;
		// make a temp for the transpose matrix
		double[][] temp = new double[A[0].length][A.length];

		int row1, row2, col1, col2;

		// copy the elements from the original matrix to the transpose matrix
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {

				temp[j][i] = A[i][j];

			}// end j for
		}// end i for

		row1 = A.length;
		row2 = temp.length;
		col1 = A[0].length;
		col2 = temp[0].length;

		return temp;

	}// end transpose

	/**
	 * Adding two matrices for calculating the covariance
	 * 
	 * @param a
	 *            the current matrix
	 * @return the sum matrix
	 */
	public double[][] addMatrixM1(double[][] a) {

		for (int c = 0; c < a.length; c++) {

			for (int d = 0; d < a[0].length; d++) {

				// add up all of the sum matrices together for the covariance
				sumMatrixM1[c][d] += a[c][d];

			}// end d for

		}// end c for

		return sumMatrixM1;

	}// end addMatrix

	/**
	 * Adding two matrices for calculating the covariance
	 * 
	 * @param a
	 *            the current matrix
	 * @return the sum matrix
	 */
	public double[][] addMatrixM2(double[][] a) {

		for (int c = 0; c < a.length; c++) {

			for (int d = 0; d < a[0].length; d++) {

				// add up all of the sum matrices together for the covariance
				sumMatrixM2[c][d] += a[c][d];

			}// end d for

		}// end c for

		return sumMatrixM2;

	}// end addMatrix

	/**
	 * check the rows and columns for multiplication
	 * 
	 * @param col1
	 *            the column of the first matrix
	 * @param row2
	 *            the row of the second matrix
	 * @return true if the column does not match, false if the do match
	 */
	public boolean checkMultiplication(int col1, int row2) {

		if (row2 != col1) {

			System.out.println("The col of the first matrix is : " + col1);
			System.out.println("The row of the second matrix is : " + row2);

			System.err
					.println("\nThe col of matrix1 does not match the row matrix2");

			return true;

		}// end if

		return false;

	}// end checkMultiplication

	/**
	 * Computes the scalar multiplication for a matrix and a number
	 * 
	 * @param a
	 *            the matrix
	 * @param b
	 *            the number to be multiplied
	 * @return the proudct matrix
	 */
	public double[][] scalarMultiplication(double[][] a, double b) {

		productMatrix = new double[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				productMatrix[i][j] = a[i][j] * b;

			}// end j for
		}// end i for

		return productMatrix;

	}// end scalarMultiplication

	/**
	 * Uses guassian to solve for the determinant
	 * 
	 * @param a
	 *            the 2D array
	 * @return the determinant matrix
	 */
	public double[][] determinant(double[][] a) {

		double[][] A = a;
		int n;
		r = 0;
		determinantNum = 0;

		n = A.length;
		// for j = 1 to n do
		for (int j = 0; j < n; j++) {

			int p = j;

			double max = Math.abs(A[p][0]);

			// find the pivot number
			// compute the pivot index j <= p <= n
			for (int i = j + 1; i < n; i++) {

				if (Math.abs(A[i][j]) > max) {
					max = A[i][j];

					p = i;

				}// end if

			}// end i for

			// check to see if there is a solution or not
			// if Cpj = 0, set E = 0 and exit.
			if (A[p][j] == 0.0) {

				determinantNum = 0;

				System.err.println("No such solution");

				System.exit(0);

			}// end if

			// switch the rows
			// if p > j, interchange rows p and j
			if (p > j) {

				double temp = 0;

				for (int i = 0; i < A[1].length; i++) {

					// copy the current element to temp
					temp = A[j][i];
					// copy the temp to another temp
					double temp2 = temp;
					// copy
					A[j][i] = A[p][i];
					double temp3 = A[p][i];
					A[p][i] = temp;

					System.out.println();

				}// end i for

				r = r + 1;

			}// end if

			// for each i != j, subtract Cij times row j from row i
			for (int i = 0; i < A.length; i++) {

				if (i > j) {

					// set the multiple of the current element
					double multiple = A[i][j];

					for (int k = 0; k < A[1].length; k++) {

						// divide the element with A[j][j] then subtract from
						// the current element
						A[i][k] = (A[i][k] - (multiple * A[j][k] / A[j][j]));

					}// end k for
				}// end if
			}// end i for

		}// end j for

		return A;
	}// end guassian

	/**
	 * Used to find the determinant
	 * 
	 * @param a
	 *            the 2D array
	 * @param b
	 *            the count
	 * @return the determinant
	 */
	public double findDelta(double[][] a, double b) {

		double determinant = 0.0, product = 1.0;

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				if (i == j) {

					// multiply the current element and put into product
					product *= a[i][j];

				}// end if

			}// end j for

		}// end i for

		determinant = Math.pow(-1, b) * product;

		return determinant;

	}// end findDelta

	/**
	 * Make the identity matrix
	 * 
	 * @param a
	 *            the covariance matrix
	 * @return the identity of the covariance matrix
	 */
	public double[][] makeIden(double[][] a) {

		idenTemp = new double[a.length][a[0].length];

		for (int k = 0; k < a.length; k++) {

			for (int l = 0; l < a[k].length; l++) {

				if (k == l) {

					// put 1s in the diagnol spots for the identity matrix
					idenTemp[k][l] = 1.0;

				}// end if

			}// end j for

		}// end i for

		return idenTemp;

	}// end makeIden

	/**
	 * Creates the augmented matrix
	 * 
	 * @param m
	 *            the first matrix
	 * @param n
	 *            the second matrix
	 * @return the augmneted matrix
	 */
	public double[][] createAugmentedMatrix(double[][] m, double[][] n) {

		constantLength = n.length;
		inverConLength = 2 * constantLength;

		// create the augmented matrix for the coefficient matrix and the
		// solution matrix

		augmentedMatrix = new double[constantLength][inverConLength];

		// create the augmented matrix
		// get the first matrix
		for (int i = 0; i < constantLength; i++) {

			for (int j = 0; j < constantLength; j++) {

				// add in the coefficient matrix to the augmented matrix

				augmentedMatrix[i][j] = m[i][j];

			}// end j for
		}// end i for

		int colCount = constantLength;

		// add the solutions to the the augmented matrix

		for (int i = 0; i < constantLength; i++) {

			for (int j = 0; j < constantLength; j++) {

				// add the constants to the end of the matrix
				augmentedMatrix[i][j + constantLength] = n[i][j];

			}// end j for

		}// end i for

		System.out.println();
		System.out.println(augmentedMatrix.length + " by "
				+ augmentedMatrix[0].length + " Matrix");

		return augmentedMatrix;

	}// end createAugmentedMatrix

	/**
	 * Computers the inverse for amtrix
	 * 
	 * @param a
	 *            the orginal matrix
	 * @return the inverted matrix
	 */
	public double[][] inverse(double[][] a) {

		double[][] C = a;
		int E, n;
		E = 1;
		n = C.length;
		// for j = 1 to n do
		for (int j = 0; j < n; j++) {

			int p = j;

			// take the absolute value of the current element and make it the
			// max
			double max = Math.abs(C[p][0]);

			// find the pivot number
			// compute the pivot index j <= p <= n
			for (int i = j + 1; i < n; i++) {

				if (Math.abs(C[i][j]) > max) {
					max = C[i][j];

					// make p = i if a new max is found
					p = i;

				}// end if

			}// end i for

			// check to see if there is a solution or not
			// if Cpj = 0, set E = 0 and exit.
			if (C[p][j] == 0.0) {

				E = 0;

				System.err.println("No such solution");

			}// end if

			// switch the rows
			// if p > j, interchange rows p and j
			if (p > j) {

				double temp = 0;

				for (int i = 0; i < C[1].length; i++) {

					// swaps the current row with the next row
					// takes the current element of the row
					temp = C[j][i];
					// stores it into a temp
					double temp2 = temp;
					C[j][i] = C[p][i];
					// copies the next pivot row's element and switches it with
					// the current element in j row
					double temp3 = C[p][i];
					// takes the current element and switch it with the pivot
					// row
					C[p][i] = temp;

				}// end i for

			}// end if

			// divide row j by the pivot
			double divide = a[j][j];

			for (int i = 0; i < C[j].length; i++) {

				// divide all of the row elements by the divider
				C[j][i] /= divide;

			}// end i for

			// for each i != j, subtract Cij times row j from row i
			for (int i = 0; i < C.length; i++) {

				if (i != j) {

					// the current element is the multiple
					double multiple = C[i][j];

					for (int k = 0; k < C[1].length; k++) {

						C[i][k] = (C[i][k] - (multiple * C[j][k]));

					}// end k for
				}// end if
			}// end i for

		}// end j for

		return C;
	}// end guassJordanSolve

	/**
	 * patitions the inverse from the inputed matrix
	 * 
	 * @param a
	 *            the inputed matrix
	 */
	public double[][] partitionInverse(double[][] a) {

		// initialize the matrices with the new length
		double[][] A = new double[a[0].length - constantLength][constantLength];
		double[][] inverse = new double[a[0].length - constantLength][constantLength];

		int k = 0;

		for (int i = 0; i < a.length; i++) {

			for (int j = a[0].length - a.length; j < a[0].length; j++) {

				// partition the inverse matrix from the identity matrix
				A[i][j - a.length] = a[i][j];

			}// end j for

		}// end i for

		return A;

	}// end print

	/**
	 * Copies the matrix of the inputed
	 * 
	 * @param a
	 *            the inputed Matrix
	 * @return the copied matrix
	 */
	public double[][] copyMatrix(double[][] a) {

		temp = new double[a.length][a[0].length];

		for (int q = 0; q < a.length; q++) {

			for (int w = 0; w < a[0].length; w++) {

				// copy the current element to the temp array
				temp[q][w] = a[q][w];

			}// end w for
		}// end q for

		return temp;

	}// end copy matrix

	/**
	 * g function
	 * 
	 * @param mean
	 * @param inverseTemp
	 * @param vectorTemp
	 * @param determinant
	 * @return the discriminant
	 */
	public double discriminant1(double[][] mean, double[][] inverseTemp,
			double[][] vectorTemp, double determinant) {

		double secondHalf = 0.0;
		double getDouble = 0.0;
		double[][] tTemp = new double[mean.length][mean[0].length];
		double[][] sTemp = new double[mean.length][mean[0].length];
		double[][] vTemp = new double[vectorTemp.length][vectorTemp[0].length];
		double[][] temp = new double[1][2];
		double[][] temp2 = new double[1][1];

		// first half
		// takes the transpose of the current vector so that the dimensions
		// equal
		// subtracts the transposed vector with the mean
		sTemp = subtraction(transpose(vectorTemp), mean);
		// multiplies the previous result with the -1/2 scalar
		tTemp = scalarMultiplication(sTemp, -0.5);
		// multiplies the previous result with the inverse of the covariance
		temp = multiply(tTemp, inverseTemp);
		// subtracts the previous result the mean that is transposed to match
		// the dimensions
		vTemp = subtraction(vectorTemp, transpose(mean));
		// multiply the previos result with
		temp2 = multiply(temp, vTemp);

		// second half
		secondHalf = ((0.5) * Math.log(determinant));

		getDouble = temp2[0][0] - secondHalf;

		return getDouble;

	}// end discriminant

	/**
	 * Multiplies two matrices
	 * 
	 * @param a
	 *            the first matrix
	 * @param b
	 *            the second matrix
	 * @return the product of the matrices
	 */
	public double[][] multiply(double[][] a, double[][] b) {

		// multiply
		// checks to see if the column of the first matrix matches the row of
		// the second matrix
		if (checkMultiplication(a[0].length, b.length) == true) {

			System.err
					.println("The column of matrix 1 does not match the row of matrix 2");

		} else

			productMatrix = new double[a.length][b[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < b[0].length; j++) {

				for (int k = 0; k < a[0].length; k++) {

					// takes the element in the first matrix and multiplies it
					// with the second matrix
					productMatrix[i][j] += a[i][k] * b[k][j];

				}// end k for
			}// end j for
		}// end i for

		return productMatrix;

	}// end multiply

	/**
	 * Subtracts two matrices
	 * 
	 * @param a
	 *            the first matrix
	 * @param b
	 *            the second matrix
	 * @return the output of the matrix
	 */
	public double[][] subtraction(double[][] a, double[][] b) {

		double[][] subSum = null;

		// runs the check to see if the dimentsion of the two matrices matches
		// each other
		if (check(a.length, a[0].length, b.length, b[0].length) == true) {

			System.err.println("The matrices do not have the same dimensions");

		} else

			subSum = new double[a.length][a[0].length];

		for (int c = 0; c < a.length; c++) {

			for (int d = 0; d < a[0].length; d++) {

				// takes the current element of the first matrix subtracts it
				// with the element in the second matrix
				subSum[c][d] = a[c][d] - b[c][d];

			}// end d for

		}// end c for

		return subSum;

	}// end subtraction

	/**
	 * check the rows and columns to see if they match
	 * 
	 * @param row1
	 *            the row of the first matrix
	 * @param col1
	 *            the column of the first matrix
	 * @param row2
	 *            the row of the second matrix
	 * @param col2
	 *            the column of the second matrix
	 * @return
	 */
	public boolean check(int row1, int col1, int row2, int col2) {

		if (row1 != row2 || col1 != col2) {

			System.err.println("\nThe dimensions do no match");

			if (row1 != row2) {

				System.err.println("\nMatrix 1 Row : " + row1
						+ " | Matrix 2 Row : " + row2);

				return true;

			}// end if
			if (col1 != col2) {

				System.err.println("\nMatrix 1 Col : " + col1
						+ " | Matrix 2 Col : " + col2);

				return true;

			}// end i

		}// end if

		// check multiplication operations

		return false;

	}// end checkl

	/**
	 * Computers the contour boundary
	 * 
	 * @param mean1
	 * @param inverse1
	 * @param determinant1
	 * @param mean2
	 * @param inverse2
	 * @param determinant2
	 * @return the contour points
	 */
	public void contour(double[][] mean1, double[][] inverse1,
			double determinant1, double[][] mean2, double[][] inverse2,
			double determinant2) {

		double output = 0.0;
		double maxX = 0.0, maxY = 0.0, minX = 9999, minY = 9999;
		double dx = 0.0, dy = 0.0;
		double x = 0.0, y = 0.0;
		double[][] temp, temp2, vector;
		double n = 50;

		// iterate through the class 1 arrays to see which x is the maxX
		// which x is the minX
		// which y is the maxY
		// which y is the minY
		for (int k = 0; k < arrays.size(); k++) {

			temp = arrays.get(k);
			temp2 = arrays2.get(k);

			if (temp[0][0] > maxX) {

				maxX = temp[0][0];

			} else if (temp[0][0] < minX) {

				minX = temp[0][0];

			}// end else if

			if (temp[1][0] > maxY) {

				maxY = temp[1][0];

			}// end if
			if (temp[1][0] < minY) {

				minY = temp[1][0];

			}// end if

			// ///////////////////////////////////////

			// iterate through the class 2 arrays to see which x is the maxX
			// which x is the minX
			// which y is the maxY
			// which y is the minY
			if (temp2[0][0] > maxX) {

				maxX = temp2[0][0];

			} else if (temp2[0][0] < minX) {

				minX = temp2[0][0];

			}// end else if

			if (temp2[1][0] > maxY) {

				maxY = temp2[1][0];

			} else if (temp2[1][0] < minY) {

				minY = temp2[1][0];

			}// end if

		}// end k for

		// calculating dx and dy
		dy = (maxY - minY) / n;
		dx = (maxX - minX) / n;

		for (int i = 0; i <= n; i++) {

			// iterating to the x point
			x = minX + (i * dx);

			for (int j = 0; j <= n; j++) {

				// iterating to the y point
				y = minY + (j * dy);
				// create a new vector with the new comouted points
				vector = new double[2][1];
				vector[0][0] = x;
				vector[1][0] = y;

				// uses both g functions to calculate the discriminant
				// takes the absolute value
				output = Math.abs(discriminant1(mean1, inverse1, vector,
						determinant1)
						- discriminant1(mean2, inverse2, vector, determinant2));

				if (output < 0.05) {

					// add the vectors are less then episilon to arrays2
					// arraylist
					arrays3.add(vector);

				}// end if

			}// end j for

		}// end i for

	}// end contour

	/**
	 * Computers the condition number for the matrix
	 * 
	 * @param a
	 *            the vector
	 * @param b
	 *            its inverse
	 * @return the condition number
	 */
	public double conditionNumber(double[][] a, double[][] b) {

		double conditionNumber = 0.0;
		double rowSum = 0.0, max1 = 0.0, max2 = 0.0;

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				// add up the elements in the row
				rowSum += Math.abs(a[i][j]);

				if (rowSum > max1) {

					// if the current rowSum is greater than max1, set max1 as
					// the rowSum
					max1 = rowSum;

				}// end if

			}// end j for

			// reset the rowSum for the next row
			rowSum = 0.0;

		}// end i for

		rowSum = 0.0;

		for (int i = 0; i < b.length; i++) {

			for (int j = 0; j < b[0].length; j++) {

				// add up the elements in the row
				rowSum += Math.abs(b[i][j]);

				if (rowSum > max2) {

					// if the current rowSum is greater than max2, set max2 as
					// the rowSum
					max2 = rowSum;

				}// end if

			}// end j for

			// reset the row sum for the next row
			rowSum = 0.0;

		}// end i for
		System.out.println("max1: " + max1);
		System.out.println("max2: " + max2);

		// multiplies the max1 with max2 to get the condition number
		conditionNumber = max1 * max2;

		return conditionNumber;

	}// end conditionNumber

	/**
	 * solve using the Guass-Jordan algorithm
	 * 
	 * @param a
	 *            the augmented matrix
	 * @return the solved matrix
	 */
	public double[][] guassJordan(double[][] a) {

		double[][] C = a;
		int E = 1;
		int n = C.length;
		// for j = 1 to n do
		for (int j = 0; j < n; j++) {

			// set p equal to the current column
			int p = j;

			// take the absolute value of the current element and make it the
			// max
			double max = Math.abs(C[p][0]);

			// find the pivot number
			// compute the pivot index j <= p <= n
			for (int i = j + 1; i < n; i++) {

				if (Math.abs(C[i][j]) > max) {
					max = C[i][j];

					// set p = i if a new max has been found in the column
					p = i;

				}// end if

			}// end i for

			// check to see if there is a solution or not
			// if Cpj = 0, set E = 0 and exit.
			if (C[p][j] == 0.0) {

				E = 0;
				// if E appears to be zero print out no solution
				System.err.println("No such solution");

			}// end if

			// switch the rows
			// if p > j, interchange rows p and j
			if (p > j) {

				double temp = 0;

				for (int i = 0; i < C[1].length; i++) {

					// swaps the current row with the next row
					// takes the current element of the row
					temp = C[j][i];
					// stores it into a temp
					double temp2 = temp;
					C[j][i] = C[p][i];
					// copies the next pivot row's element and switches it with
					// the current element in j row
					double temp3 = C[p][i];
					// takes the current element and switch it with the pivot
					// row
					C[p][i] = temp;

				}// end i for

			}// end if

			// divide row j by the pivot
			double divide = a[j][j];

			for (int i = 0; i < C[j].length; i++) {

				// divides the row by the current element where j = j
				C[j][i] /= divide;

			}// end i for

			// for each i != j, subtract Cij times row j from row i
			for (int i = 0; i < C.length; i++) {

				if (i != j) {

					double multiple = C[i][j];

					for (int k = 0; k < C[1].length; k++) {

						// multiplies the current element in j row with the
						// multiple
						// then subtracts it from the current row's element
						C[i][k] = (C[i][k] - (multiple * C[j][k]));

					}// end k for
				}// end if
			}// end i for

		}// end j for

		return C;
	}// end guassJordanSolve

	/**
	 * Prints the output of the guassJordan algorithm
	 * 
	 * @param a
	 */
	public void printOutput(double[][] a) {

		// print the solution to the elimination

		System.out.println("Final Output: ");

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				System.out.print(a[i][j] + " ");

			}// end j for

			System.out.println();
		}// end i for

		System.out.println("\nThe solution is : ");

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				if (j == a.length) {

					// prints out the solutions
					System.out.println(a[i][a.length]);

				}// end if
			}// end j for
		}// end i for
	}// end print output

	/**
	 * Partitions the coefficient from the matrix A
	 * 
	 * @param a
	 *            the matrix A
	 * @return the coefficient matrix
	 */
	public double[][] partitionSystems(double[][] a) {

		double[][] partitionSys = new double[a.length][a[0].length - 1];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length - 1; j++) {

				// takes out the coefficients from the constants of matrix A
				partitionSys[i][j] = a[i][j];

			}// end j for

		}// end i for

		return partitionSys;

	}// end partitionSystems

}// end class matrix