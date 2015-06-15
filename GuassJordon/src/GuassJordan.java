import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.StringTokenizer;

public class GuassJordan {

	public int choice, constantLength, inverConLength, countTransposeRow,
			countTransposeCol;
	public double scalar, determinantNum, r;
	public double[][] inverseM1, inverseM2, idenMatrixM1, idenMatrixM2,
			augmentedMatrix, productMatrix, sumMatrixM1, sumMatrixM2,
			determinantM1, determinantM2, temp;
	public double[][] vector, vector2, tempVector, transposeMatrix,
			covarianceMatrixM1, covarianceMatrixM2, identityMatrix,
			augMatrixM1, augMatrixM2, idenTemp;
	public boolean matrix1Created;
	public Scanner input2;
	public ArrayList<double[][]> arrays = new ArrayList<double[][]>(),
			arrays2 = new ArrayList<double[][]>();

	public double[][] test;
	public ArrayList<Integer> testInts;

	public static void main(String[] args) {

		GuassJordan obj = new GuassJordan();
		obj.run();

	}// end main

	public void run() {

		double[][] temp = createMatrix();
		sumMatrixM1 = new double[2][2];
		sumMatrixM2 = new double[2][2];
		covarianceMatrixM1 = new double[2][2];
		covarianceMatrixM2 = new double[2][2];
		determinantM1 = new double[2][2];
		determinantM1 = new double[2][2];
		double[][] tempCovM1 = new double[2][2];
		double[][] tempCovM2 = new double[2][2];

		computeMean(temp);

		// calculate the covariance of M1
		for (int i = 0; i < arrays.size(); i++) {

			tempVector = arrays.get(i);
			transposeMatrix = transpose(tempVector);
			addMatrixM1(transposeMatrix);

		}// end i for

		// calculate the covariance of M2
		for (int i = 0; i < arrays2.size(); i++) {

			tempVector = arrays2.get(i);
			transposeMatrix = transpose(tempVector);
			addMatrixM2(transposeMatrix);

		}// end i for

		// multiply the sum matrix M1 by the scalar
		covarianceMatrixM1 = scalarMultiplication(sumMatrixM1, scalar);
		// multiply the sum matrix M2 by the scalar
		covarianceMatrixM2 = scalarMultiplication(sumMatrixM2, scalar);

		tempCovM1 = copyMatrix(covarianceMatrixM1);
		tempCovM2 = copyMatrix(covarianceMatrixM2);

		// solve for the determinant using gaussian implementation
		determinantM1 = determinant(covarianceMatrixM1);
		determinantM2 = determinant(covarianceMatrixM2);

		System.out.println("\nDeterminant of M1 = "
				+ findDelta(determinantM1, r));
		System.out.println("\nDeterminant of M2 = "
				+ findDelta(determinantM2, r));

		idenMatrixM1 = makeIden(tempCovM1);
		idenMatrixM2 = makeIden(tempCovM2);

		augMatrixM1 = createAugmentedMatrix(tempCovM1, idenMatrixM1);
		augMatrixM2 = createAugmentedMatrix(tempCovM2, idenMatrixM2);

		inverseM1 = inverse(augMatrixM1);
		inverseM2 = inverse(augMatrixM2);

		System.out.println("\nPrinting the inverse matrices: ");
		System.out.println("M1: ");
		printMatrix(inverseM1);
		System.out.println("M2: ");
		printMatrix(inverseM2);
		System.out.println("\n...The inverse has been computed");

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

		Scanner sc = new Scanner(System.in);
		String fileName = sc.next();

		double number = 0;
		int row = 0, col = 0;

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

				col = st.countTokens();
				row++;

			}// end while

			createdMatrix = new double[row][col];

			System.out.println("Created a matrix of size : " + row + " by "
					+ col);

			fin = new FileReader(file);
			matrixFile3 = new BufferedReader(fin);

			while ((line = matrixFile3.readLine()) != null) {

				StringTokenizer st = new StringTokenizer(line);

				for (j = 0; j < col; j++) {

					if (i != 0) {

						number = Double.parseDouble(st.nextToken());
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
			for (double[][] c : arrays) {

				for (int p = 0; p < c.length; p++) {

					for (int q = 0; q < c[0].length; q++) {

						System.out.print("[" + c[p][q] + "]");

					}// end q
				}// end p

				System.out.println();
			}// end for

			System.out.println("Array 1 size: " + arrays.size());
			System.out.println("Array 2 size: " + arrays2.size());
		} catch (IOException e) {
			System.err.println(e);
		}// end catch

		return temp;

	}// end createMatrix

	/**
	 * prints out the matrices
	 * 
	 * @param a
	 *            the matrix to print
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

	public void debug(double[][] d) {

		System.out.println();
		System.out.println("Debugging...");

		for (int i = 0; i < d.length; i++) {

			for (int k = 0; k < d[1].length; k++) {

				System.out.print("[" + d[i][k] + "]");
			}// end k for

			System.out.println();
		}// end i for
	}// end debug

	public void split(double[][] a) {

		System.out.println("printing...");
		System.out.println("The length of a is : " + a[0].length);

		double[][] A = new double[a[0].length - constantLength][constantLength];
		double[][] inverse = new double[a[0].length - constantLength][constantLength];

		System.out.println("A: " + A.length + " by " + A[0].length + " Matrix");
		System.out.println("Inverse: " + inverse.length + " by "
				+ inverse[0].length + " Matrix");

		System.out.println("a: " + a.length + " by " + a[0].length + " Matrix");

		System.out.println("\nThe ivnerse matrix is : ");

		int k = 0;

		for (int i = 0; i < a.length; i++) {

			for (int j = a[0].length - a.length; j < a[0].length; j++) {

				A[i][j - a.length] = a[i][j];

				System.out.print("[" + A[i][j - a.length] + "]");

			}// end j for

			System.out.println();

		}// end i for

		System.out.println();
		System.out.println(a.length + " by " + a[0].length + " Matrix");

	}// end split

	public void computeMean(double[][] a) {

		double[][] A = a;
		double[][] mean1, mean2;
		double sumX1 = 0, sumY1 = 0, sumX2 = 0, sumY2 = 0, multiple;
		scalar = Math.pow(A.length, -1);
		System.out.println("Scalar: " + scalar);

		mean1 = new double[1][A[0].length / 2];
		mean2 = new double[1][A[0].length / 2];

		// mean for x1

		for (int i = 0; i < A.length; i++) {

			for (int j = 0; j < A[0].length; j++) {

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

		mean1[0][0] = sumX1 * scalar;
		mean1[0][1] = sumX2 * scalar;

		System.out.println("\nMean of m1: ");

		for (int k = 0; k < mean1.length; k++) {

			for (int l = 0; l < mean1[0].length; l++) {

				System.out.print("[" + mean1[k][l] + "]");

			}// end l for

		}// end k for

		if (A[0].length >= 4) {
			mean2[0][0] = sumY1 * scalar;
			mean2[0][1] = sumY2 * scalar;

			System.out.println("\nMean of m2: ");

			for (int k = 0; k < mean2.length; k++) {

				for (int l = 0; l < mean2[0].length; l++) {

					System.out.print("[" + mean2[k][l] + "]");

				}// end l for

			}// end k for
		}// end if
	}// end compute mean

	public double[][] transpose(double[][] a) {

		double[][] A = a;
		// each 2x1 matrix
		double[][] temp = new double[A[0].length][A.length];

		int row1, row2, col1, col2;

		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {

				temp[j][i] = A[i][j];

			}// end j for
		}// end i for

		row1 = A.length;
		row2 = temp.length;
		col1 = A[0].length;
		col2 = temp[0].length;

		// multiply transpose
		if (checkMultiplication(col1, row2) == true) {

			System.err.println("Column of the first matrix: " + col1
					+ "does not match with the row of the first matrix: "
					+ row2);

		} else

			productMatrix = new double[row1][col2];

		for (int i = 0; i < row1; i++) {

			for (int j = 0; j < col2; j++) {

				for (int k = 0; k < col1; k++) {
					productMatrix[i][j] += A[i][k] * temp[k][j];

				}// end k for
			}// end j for
		}// end i for

		// System.out.println("Solution to transpose");

		return productMatrix;

	}// end transpose

	public double[][] addMatrixM1(double[][] a) {

		for (int c = 0; c < a.length; c++) {

			for (int d = 0; d < a[0].length; d++) {

				sumMatrixM1[c][d] += a[c][d];

			}// end d for

		}// end c for

		return sumMatrixM1;

	}// end addMatrix

	public double[][] addMatrixM2(double[][] a) {

		for (int c = 0; c < a.length; c++) {

			for (int d = 0; d < a[0].length; d++) {

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

		r = 0;
		int n;
		determinantNum = 0;
		r = 0;
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

					System.out.println("\nCompute the pivot");

					System.out.print("\nThe max is : " + max);

					p = i;

					System.out.println("\nThe p pivot index is now: " + p);

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

				System.out.println("\nSwap the rows");
				System.out.println("if " + p + " is greater then " + j);

				double temp = 0;

				for (int i = 0; i < A[1].length; i++) {

					temp = A[j][i];
					double temp2 = temp;
					A[j][i] = A[p][i];
					double temp3 = A[p][i];
					A[p][i] = temp;

					System.out.println(temp2 + " switch with " + temp3);

					System.out.println();

				}// end i for

				r = r + 1;
				System.out.println("The total r is : " + r);

			}// end if

			// for each i != j, subtract Cij times row j from row i
			for (int i = 0; i < A.length; i++) {

				if (i > j) {

					double multiple = A[i][j];

					System.out
							.println("For each i != j, subtract Cij times row j from row i");
					System.out.println("The i : " + j);
					System.out.println("The j : " + i);

					System.out.println("The multiple is : " + A[i][j]);

					for (int k = 0; k < A[1].length; k++) {

						System.out.println("Multiply then subract : " + A[i][k]
								+ " - " + multiple + " * " + A[j][k]);

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

		System.out.println("Finding the determinant");

		double determinant = 0.0, product = 1.0;

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				if (i == j) {

					product *= a[i][j];
					System.out.println(product);

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

					idenTemp[k][l] = 1.0;

				}// end if

			}// end j for

		}// end i for

		return idenTemp;

	}// end makeIden

	/**
	 * 
	 * @param m
	 * @param n
	 * @return
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

				System.out.println("The i is : " + i);
				System.out.println("The j is : " + j);

				augmentedMatrix[i][j + constantLength] = n[i][j];

			}// end j for

		}// end i for

		// print the augmented matrix

		System.out.println("Displaying Augmented Matrix...");

		if (augmentedMatrix.length != 0) {

			System.out.println("_____Augmented Matrix_____");

			for (int k = 0; k < augmentedMatrix.length; k++) {

				for (int j = 0; j < augmentedMatrix[0].length; j++) {

					System.out.print(augmentedMatrix[k][j] + "\t");

				}// end j for

				System.out.println();

			}// end i for

			System.out.println();
			System.out.println(augmentedMatrix.length + " by "
					+ augmentedMatrix[0].length + " Matrix");

		}// end if

		return augmentedMatrix;

	}// end createAugmentedMatrix

	public double[][] inverse(double[][] a) {

		double[][] C = a;
		int E, n;
		E = 1;
		n = C.length;
		// for j = 1 to n do
		for (int j = 0; j < n; j++) {

			int p = j;

			double max = Math.abs(C[p][0]);

			// find the pivot number
			// compute the pivot index j <= p <= n
			for (int i = j + 1; i < n; i++) {

				if (Math.abs(C[i][j]) > max) {
					max = C[i][j];

					System.out.println("\nCompute the pivot");

					System.out.print("\nThe max is : " + max);

					p = i;

					System.out.println("\nThe p pivot index is " + p);

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

				System.out.println("\nSwap the rows");
				System.out.println("if " + p + " is greater then " + j);

				double temp = 0;

				for (int i = 0; i < C[1].length; i++) {

					temp = C[j][i];
					double temp2 = temp;
					C[j][i] = C[p][i];
					double temp3 = C[p][i];
					C[p][i] = temp;

					System.out.println(temp2 + " switch with " + temp3);

					System.out.println();

				}// end i for

			}// end if

			// divide row j by the pivot
			double divide = a[j][j];

			System.out.println("The divisor is : " + divide);

			for (int i = 0; i < C[j].length; i++) {

				System.out.println("Divide : " + C[j][i] + " / " + divide);

				C[j][i] /= divide;

			}// end i for

			// for each i != j, subtract Cij times row j from row i
			for (int i = 0; i < C.length; i++) {

				if (i != j) {

					double multiple = C[i][j];

					System.out
							.println("For each i != j, subtract Cij times row j from row i");

					System.out.println("The multiple is : " + C[i][j]);

					for (int k = 0; k < C[1].length; k++) {

						System.out.println("Multiply then subract : " + C[i][k]
								+ " - " + multiple + " * " + C[j][k]);

						C[i][k] = (C[i][k] - (multiple * C[j][k]));

					}// end k for
				}// end if
			}// end i for

		}// end j for

		return C;
	}// end guassJordanSolve

	/**
	 * prints out the matrices
	 * 
	 * @param a
	 *            the matrix to print
	 */
	public void printInverse(double[][] a) {

		System.out.println("printing...");
		System.out.println("The length of a is : " + a[0].length);

		double[][] A = new double[a[0].length - constantLength][constantLength];
		double[][] inverse = new double[a[0].length - constantLength][constantLength];

		System.out.println("A: " + A.length + " by " + A[0].length + " Matrix");
		System.out.println("Inverse: " + inverse.length + " by "
				+ inverse[0].length + " Matrix");

		System.out.println("a: " + a.length + " by " + a[0].length + " Matrix");

		System.out.println("\nThe ivnerse matrix is : ");

		int k = 0;

		for (int i = 0; i < a.length; i++) {

			for (int j = a[0].length - a.length; j < a[0].length; j++) {

				A[i][j - a.length] = a[i][j];

				System.out.print("[" + A[i][j - a.length] + "]");

			}// end j for

			System.out.println();

		}// end i for

		System.out.println();
		System.out.println(a.length + " by " + a[0].length + " Matrix");

		// System.out.println("The constant length is: " + constantLength);
		//
		// System.out.println("\nA^-1 = ");
		//
		// for (int i = 0; i < a.length; i++) {
		//
		// for (int j = 0; j < a[0].length; j++) {
		//
		// System.out.print(a[i][j + constantLength]);
		//
		// }// end j for
		//
		// System.out.println();
		//
		// }// end i for

	}// end print

	/**
	 * Copies the matrix of the inputted
	 * 
	 * @param a
	 *            the inputted Matrix
	 * @return the copied matrix
	 */
	public double[][] copyMatrix(double[][] a) {

		temp = new double[2][2];

		for (int q = 0; q < a.length; q++) {

			for (int w = 0; w < a[0].length; w++) {

				temp[q][w] = a[q][w];

			}// end w for
		}// end q for

		return temp;

	}// end copy matrix

}// end class matrix