import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Scanner;
import java.util.StringTokenizer;

/**
 * 
 */

/**
 * Class to test # number 1 and 2 of the test
 * 
 * @author Kavy Rattana
 * 
 */
public class test2part1 {

	public int choice, r;
	double determinantNum;
	public String fileName;
	public Vector[] arrayPoints;
	Vector meanV;
	public double[][] covariance, productMatrix, covarianceMatrix;
	public double determinant, scalar;
	public double[][] mean, temp, sumMatrix;

	/**************************************************************************
	 * Running the program The text file should be called "data2.txt"
	 * 
	 **************************************************************************/

	public static void main(String[] args) throws IOException {

		test2part1 obj = new test2part1();
		obj.run();

	}// end main

	public void run() throws IOException {

		// initialize the variables
		covarianceMatrix = new double[2][2];
		determinant = 0.0;

		createMatrix();

		printMatrix(arrayPoints);

		System.out.println("Choose number to run:");
		System.out.println("1. Compute Mean, Covariance, Eigenvalues");
		System.out.println("2. Leveriier, Power Method");

		Scanner sc = new Scanner(System.in);
		choice = sc.nextInt();

		// Problem 1. a
		// find the mean
		if (fileName.equals("data.txt") && choice == 1) {

			System.out.println("\nThe Mean:");
			double[][] meanMatrix = computeMean(arrayPoints);
			printMatrix2D(mean);

			meanV = new Vector(mean[0][0], mean[0][1]);

			// initializing variables
			double[][] covMatrix = new double[2][1];
			double[][] transposeMatrix = new double[2][1];
			double[][] sumMatrix = new double[2][2];
			double[][] tempTrans = new double[2][1];
			double[][] tempE = new double[2][2];
			double[][] eigenValue = new double[2][2];
			double[][] determinant = new double[2][2];
			double[] eigen = new double[2];
			double[][] identity = new double[2][2];

			// compute the covariance

			for (int i = 0; i < arrayPoints.length; i++) {

				// subtract the mean from each array point
				covMatrix[0][0] = (arrayPoints[i].getX() - meanV.getX());
				covMatrix[1][0] = (arrayPoints[i].getY() - meanV.getY());

				// compute the transpose
				tempTrans = copyMatrix(covMatrix);
				transposeMatrix = transpose(tempTrans);

				// multiply the covariance matrix with the transpose
				tempE = multiply(covMatrix, transposeMatrix);

				// add the sum matrix to temp
				sumMatrix = sumMatrix(sumMatrix, tempE);

			}// end for

			// divide each element in sum matrix by the number of points
			for (int g = 0; g < sumMatrix.length; g++) {

				for (int j = 0; j < sumMatrix[0].length; j++) {

					sumMatrix[g][j] = sumMatrix[g][j] / arrayPoints.length;

				}// end j for

				System.out.println("");

			}// end for

			// compute the trance of the sum matric
			double TRACE = trace(sumMatrix);
			System.out.println("The trace of the covariance is " + TRACE);

			System.out.println("Compute Eigen Values");

			// print out the sum matrix
			debug(sumMatrix);

			// copy the sum matrices
			determinant = copyMatrix(sumMatrix);
			eigenValue = copyMatrix(sumMatrix);

			determinant = determinant(determinant);

			
			System.out.println("Determinant");
			debug(determinant);

			// compute the eigen values
			eigenvalues(sumMatrix, identity, eigen);
			System.out.println("THE EIGEN VALUES______");
			for (int i = 0; i < eigen.length; i++) {

				System.out.println(eigen[i]);
			}

			// Computing the eigen values
			double set = 1;

			// finding the
			double LAB1 = ((eigen[0] - eigenValue[0][0]) / eigenValue[0][1]);
			double LAB2 = ((eigen[1] - eigenValue[0][0]) / eigenValue[0][1]);

			// make a new vector from based on values from LAB 1 an LAB 2
			Vector EIG1 = new Vector(set, LAB1);
			Vector EIG2 = new Vector(set, LAB2);

			// Print the eigen values
			System.out.println(EIG1.tostring() + " :: " + EIG2.tostring());

			// find unit legnths of the vectors
			Vector unitLength1 = unitv(EIG1);
			Vector unitLength2 = unitv(EIG2);

			// print the unit lengths
			System.out.println(unitLength1.tostring() + " :: "
					+ unitLength2.tostring());

		}// end if for question 1

		// Running leverier's and power method
		if (fileName.equals("data.txt") && choice == 2) {

			double[][] A1 = { { 4, 17, -24, -36 }, { 1, 0, 0, 0 },
					{ 0, 1, 0, 0 }, { 0, 0, 1, 0 } };
			double[][] B1 = { { 1, 1, 1, 1 } };

			// check the array A
			debug(A1);

			// compute Leverier method
			double[] leverierMatrix = leverier(A1);

			System.out.println("Leverier Matrix");

			for (int i = 0; i < leverierMatrix.length; i++) {

				System.out.println(" " + leverierMatrix[i]);

			}

			// compute power method
			double powerNumber = powerMethod(A1);
			System.out.println("The power function returned mu" + powerNumber);

		}// end if

	}// end run

	/**************************************************************************
	 * Creating the matrix
	 * 
	 **************************************************************************/

	/**
	 * create the matrices
	 * 
	 * @return double[][] the specific matrix
	 */
	public double[][] createMatrix() {

		double[][] vector;
		Vector point;
		int i = 0, j = 0, index = 0;
		char label = 0;
		double x = 0, y = 0;
		int row = 0, col = 0;
		double[][] createdMatrix = null;
		double temp[][] = null;

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

			arrayPoints = new Vector[tempMatrix.length];

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

					}// end if

					point = new Vector(x, y);
					arrayPoints[i - 1] = point;// array of vectors

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
	 * Compute Mean
	 * 
	 **************************************************************************/

	/**
	 * Computes the mean of the matrix
	 * 
	 * @param a
	 *            the array of vectors
	 * @return the mean vector
	 */
	public double[][] computeMean(Vector[] a) {

		double sumX = 0.0;
		double sumY = 0.0;
		double totalPoints = 0.0;

		mean = new double[1][2];

		for (int i = 0; i < a.length; i++) {

			sumX += a[i].getX();
			sumY += a[i].getY();

			System.out.println("sumX: " + sumX);
			System.out.println("sumY: " + sumY);

			totalPoints++;

		}// end i for

		mean[0][0] = sumX / totalPoints;
		mean[0][1] = sumY / totalPoints;

		return mean;

	}// end computeMean

	/**
	 * /************************************************************************
	 * ** Eigenvalues
	 * 
	 **************************************************************************/

	/**
	 * 
	 * Finds the eigenvalues of a matrix using the Jacobi Iterative method.
	 * 
	 * @param covarianceTemp
	 *            the covariance matrix
	 * @param idenMatrix
	 *            the identity matrix
	 * @param EigenTemp
	 *            the eigen matrix
	 */
	public void eigenvalues(final double covarianceTemp[][],
			double idenMatrix[][], double EigenTemp[]) {

		int n = covarianceTemp.length;
		double ATemp[][] = new double[n][n];
		double cTemp[] = new double[1];
		double sTemp[] = new double[1];

		if (covarianceTemp[0].length != n || idenMatrix.length != n
				|| idenMatrix[0].length != n || EigenTemp.length != n) {

			System.err.println("The Matrices do not match");

		}

		// Initialize
		cTemp[0] = 1.0;
		sTemp[0] = 0.0;

		for (int i = 0; i < n; i++) {

			for (int j = 0; j < n; j++) {

				idenMatrix[i][j] = 0.0;

			}// end j for

			idenMatrix[i][i] = 1.0;

		}// end i for
			// make a copy of the covariance matric to use
		ATemp = copyMatrix(covarianceTemp);

		for (int k = 0; k < n; k++) {
			for (int i = 0; i < n - 1; i++) {
				for (int j = i + 1; j < n; j++) {

					computeT(ATemp, i, j, cTemp, sTemp);
					multTrans(i, j, cTemp, sTemp, ATemp, idenMatrix);

				}// end j for
			} // end i for
		}// end k for

		// re copy the eigen values
		for (int i = 0; i < n; i++)

			EigenTemp[i] = ATemp[i][i];

	} // end eigenvalues

	/**
	 * Compute theta
	 * 
	 * @param A
	 * @param p
	 * @param q
	 * @param c
	 * @param s
	 */
	public void computeT(final double A[][], final int p, final int q,
			double c[], double s[]) {

		double theta;
		double t;

		if (A[p][q] != 0.0) {

			theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
			if (theta >= 0.0) {

				t = 1.0 / (theta + Math.sqrt(1.0 + theta * theta));

			} else {

				t = -1.0 / ((-theta) + Math.sqrt(1.0 + theta * theta));
			}

			c[0] = 1.0 / Math.sqrt(1.0 + t * t);
			s[0] = t * c[0];

		} else {

			c[0] = 1.0;
			s[0] = 0.0;
		}

	} // end computeT

	/**
	 * Multiplies three matrices together
	 * 
	 * @param A
	 *            the first matrix
	 * @param B
	 *            the second matrix
	 * @param C
	 *            the third matrix
	 */
	public void multiply(double A[][], final double B[][], double C[][]) {

		int row1 = A.length; // the row of matrix 1
		int col1 = A[0].length; // the col of matrix 1
		int row2 = B[0].length; // the col of matrix 2

		if (B.length != col1 || C.length != row1 || C[0].length != row2) {

			System.out.println("Error in Matrix.multiply, incompatible sizes");
		}

		// multiply the three matrices
		for (int i = 0; i < row1; i++)

			for (int j = 0; j < row2; j++) {

				C[i][j] = 0.0;

				for (int k = 0; k < col1; k++) {

					C[i][j] = C[i][j] + A[i][k] * B[k][j];

				}// end k for
			}// end j for

	} // end multiply

	/**
	 * Transposes several matrcies
	 * 
	 * @param p
	 * @param q
	 * @param c
	 * @param s
	 * @param A
	 * @param V
	 */
	public void multTrans(int p, final int q, double c[], double s[],
			double A[][], double V[][]) {

		int n = A.length;
		double B[][] = new double[n][n];
		double J[][] = new double[n][n];
		// reset j
		for (int i = 0; i < n; i++) {

			for (int j = 0; j < n; j++) {

				J[i][j] = 0.0;

			}// end j for

			J[i][i] = 1.0;

		}// end i for

		/* J transpose */
		J[p][p] = c[0];
		J[p][q] = -s[0];
		J[q][q] = c[0];
		J[q][p] = s[0];

		B = multiply(J, A);
		J[p][q] = s[0];
		J[q][p] = -s[0];
		multiply(B, J, A);
		multiply(V, J, B);
		V = copyMatrix(B);

	} // multTrans

	/**
	 * Finds the unit length of a vector
	 * 
	 * @param A
	 *            - the current vector
	 * @return the unit length vector
	 */
	public Vector unitv(Vector A) {

		double x = A.getX();
		double y = A.getY();
		double z = Math.sqrt(Math.pow(A.getX(), 2) + Math.pow(A.getY(), 2));

		Vector b = new Vector((x / z), (y / z));
		return b;

	}

	/**
	 * /************************************************************************
	 * ** Leverier Method
	 * 
	 **************************************************************************/

	/**
	 * Implement Leverier method to characteristics
	 * 
	 * @param A
	 *            - the matrix
	 * @return
	 */
	public double[] leverier(double[][] A) {

		// recopy the matrix to be used
		double[][] B = copyMatrix(A);
		int n = A.length;

		double[] a = new double[A.length];
		System.out.println("A is " + a.length);
		a[n - 1] = (0 - (trace(A)));

		for (int k = (n - 1); k >= 1; k--) {

			System.out.println("Trace Addition + Multiply");
			B = traceAddition(B, a[k]);
			B = multiply(A, B);

			a[k - 1] = (0 - trace(B) / ((n - k) + 1));

			System.out.println("Leverier Method: " + a[k - 1]);
		}

		return a;

	}// end leverier

	/**
	 * Adds a double to the trace of a matrix (used for leverier's method)
	 * 
	 * @param A
	 *            - the matrix being added to
	 * @param x
	 *            - the double being added
	 * @return
	 */
	public double[][] traceAddition(double[][] A, double x) {

		for (int i = 0; i < A.length; i++) {

			A[i][i] = (A[i][i] + x);

		}

		return A;

	}// end trace add

	/**
	 * /************************************************************************
	 * ** Power Method
	 * 
	 **************************************************************************/

	/**
	 * Implements the power method to find the most dominant eigenvalue
	 * 
	 * @param A
	 *            the matrix
	 * @return the value of x
	 */
	public double powerMethod(double[][] A) {

		double episilon = .00001;
		double m = (((A.length + A[0].length) * 100) + 1);
		double[][] y = new double[A.length][1];

		for (int i = 0; i < y.length; i++) {

			y[i][0] = (Math.random() * 1) + 1;

		}// end for

		int k = 0;

		double[][] x = multiply(A, y);

		// print out the value of x
		System.out.println("x: ");
		debug(x);

		double r = 1;
		double mu = 0;

		while (k < m && r > episilon) {

			double TM = 1 / magnitude(x);
			y = scalarMultiplication(x, TM);
			x = multiply(A, y);

			double[][] yt = transpose(y);
			double[][] numerator = multiply(yt, x);
			double NUMER = makeD(numerator);

			double[][] denom = multiply(yt, y);
			double DENOM = makeD(denom);

			mu = NUMER / DENOM;

			double[][] z = scalarMultiplication(y, mu);
			double[][] q = subtraction(z, x);
			r = magnitude(q);
			k++;

		}// end while

		return mu;

	}// end power method

	/**
	 * Finds the magnitude of the matrix
	 * 
	 * @param A
	 *            the matrix
	 * @return the magnitude
	 */
	public double magnitude(double[][] A) {

		double mag = 0;

		for (int i = 0; i < A.length; i++) {

			for (int j = 0; j < A[0].length; j++) {

				mag = mag + Math.pow(A[i][j], 2);

			}// end j for
		}// end i for

		mag = Math.sqrt(mag);
		return mag;

	}// end magnitude

	/**
	 * /************************************************************************
	 * ** Operations
	 * 
	 **************************************************************************/

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
	 * Turns a matrix in to double
	 * 
	 * @param a
	 *            - the 1x1 matrix
	 * @return b (a double)
	 */
	public static double makeD(double[][] a) {
		double b = 0;
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {

				b = a[i][j];
				// System.out.println("\n B found at "+b+"\n");

			}// end i for

		}// end j for
			// System.out.println("\n B returned "+b+"\n");
		return b;

	}// end find

	/**
	 * Prints the contents of a amatric
	 * 
	 * @param a
	 *            the matrix
	 */
	public void printMatrix(Vector[] a) {

		for (int i = 0; i < a.length; i++) {

			System.out.print(a[i].getX() + " " + a[i].getY());
			System.out.println();

		}// end i for

		System.out.println("Length of array: " + a.length);

	}// end printMatrix

	/**
	 * Prints the 2D matrix
	 * 
	 * @param a
	 *            the 2D matrix
	 */
	public void printMatrix2D(double[][] a) {

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				System.out.print(a[i][j] + " ");

			}// end j for

			System.out.println();

		}// end i for

	}// end printMatrix2D

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
	 * Adding two matrices for calculating the covariance
	 * 
	 * @param a
	 *            the current matrix
	 * @return the sumMatrix matrix
	 */
	public double[][] sumMatrix(double[][] a, double[][] b) {

		sumMatrix = new double[a.length][a[0].length];

		for (int c = 0; c < a.length; c++) {

			for (int d = 0; d < a[0].length; d++) {

				// add up all of the sumMatrix matrices together for the
				// covariance
				sumMatrix[c][d] = a[c][d] + b[c][d];

			}// end d for

		}// end c for

		return sumMatrix;

	}// end addMatrix

	/**
	 * 
	 * Finds the trace of the matrix
	 * 
	 * @param a
	 *            - the matrix
	 * 
	 * @return total- the trace
	 */
	public static double trace(double[][] a) {
		double total = 0;
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {
				if (i == j) {
					total = total + a[i][j];
				}// end if
			}

		}

		return total;

	}// end debug

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

}// end class
