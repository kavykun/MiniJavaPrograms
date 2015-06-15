import java.util.ArrayList;

/**
 * First part of the test
 * 
 * @author Kavy Rattana
 * 
 */
public class Test3part1 {

	double constantLength, inverConLength;
	double[][] polyMatrix, augmentedMatrix, y, idenMatrix, inverse, temp,
			productMatrix;
	ArrayList<Vector> numbers;
	double r;

	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

		Test3part1 obj = new Test3part1();
		obj.run();

	}// end main

	/**
	 * Method to run the program
	 */
	public void run() {

		numbers = new ArrayList<Vector>();
		Vector v1 = new Vector(-2, -4);
		Vector v2 = new Vector(-1, 1);
		Vector v3 = new Vector(0, -6);
		Vector v4 = new Vector(1, -1);

		numbers.add(v1);
		numbers.add(v2);
		numbers.add(v3);
		numbers.add(v4);

		polyMatrix = new double[numbers.size()][numbers.size()];
		y = new double[4][1];

		// fill in the matrix
		for (int i = 0; i < numbers.size(); i++) {

			for (int j = 0; j < numbers.size(); j++) {

				polyMatrix[i][j] = numbers.get(i).getX();

			}// end j for

		}// end i for

		for (int i = 0; i < y.length; i++) {

			for (int j = 0; j < 1; j++) {

				y[i][j] = numbers.get(i).getY();

			}// end j for

		}// end i for

		polyMatrix = polynomial(polyMatrix);

		idenMatrix = makeIden(polyMatrix);

		augmentedMatrix = createAugmentedMatrix(polyMatrix, idenMatrix);

		inverse = partitionInverse(inverse(augmentedMatrix));

		productMatrix = multiply(inverse, y);

		System.out.println("Product Matrix: ");

		debug(productMatrix);

		// debug(productMatrix);

		// Root Finding

		double[][] output;

		System.out.println("Poly Matrix ");
		debug(productMatrix);

		output = bisectionMethod(-2, -1, 0, 1, productMatrix);

		for (int i = 0; i < output.length; i++) {

			for (int j = 0; j < output[0].length; j++) {

				System.out.print(output[i][j] + "\t");

			}// end j for

			System.out.println();

		}// end i for

	}// end run

	/**
	 * Runs the bisection method
	 * 
	 * @param n1
	 * @param n2
	 * @param n3
	 * @param n4
	 * @return
	 */
	public double[][] bisectionMethod(double n1, double n2, double n3,
			double n4, double[][] matrix) {

		double[][] matrixOutput = new double[4][1];
		double[] xValues = { n1, n2, n3, n4 };
		double x0, x1, x2;
		double f0, f2;
		double eps = .000000000001;
		int x = 0;
		while (x < 3) {
			x0 = xValues[x];
			x1 = xValues[x + 1];
			for (int m = 50000, k = 1; m > 0; m--) {
				f0 = functionBisection(x0, matrix);
				do {
					x2 = (x0 + x1) / 2;

					f2 = functionBisection(x2, matrix);

					if ((f0 * f2) < 0) {
						x1 = x2;
					} else {
						x0 = x2;
						f0 = f2;
					}
					k++;
				} while (Math.abs(f2) > eps && k <= m);

				matrixOutput[x][0] = x2;
			}
			x++;

		}
		return matrixOutput;
	}

	/**
	 * This method does the lagrange polynomial with an x value given. The
	 * coefficient values are taken from the calling coefficient matrix.
	 * 
	 * @param x
	 *            the x value for polynomial
	 * @return the result
	 */
	public double functionBisection(double x, double[][] matrix) {

		System.out.println("Output");
		debug(matrix);

		double f = matrix[0][0] * Math.pow(x, 3) + matrix[1][0]
				* Math.pow(x, 2) + matrix[2][0] * x + matrix[3][0];

		return f;
	}

	// /**
	// * Find the roots using the muller method
	// */
	// public void muller() {
	//
	// double x0, x1, m, f, dm, x, x2, k, f0, f1, eps, y;
	// double h0, h1, d, c0, ca, cb, c1, c2, x3, f2, c;
	//
	// // An acceptable level of precision e>0
	// // An iteration limit m>0
	// // Interval endpoints x0 and x1 such that f(x0)*f(x1)<0
	// // Set f0 = f(x0), f1 = f(x1)
	// // Set k=1
	//
	// k = 1;
	// y = eps + 1;
	//
	// do {
	// h0 = x0 - x2;
	// h1 = x1 - x2;
	// d = h0 * h1 * (h1 - h0);
	// c0 = f2;
	//
	// ca = h1 * h1 * (f0 - f2);
	// cb = h0 * h0 * (f1 - f2);
	// c1 = (ca - cb) / d;
	//
	// ca = h0 * (f1 - f2);
	// cb = h1 * (f0 - f2);
	// c2 = (ca - cb) / d;
	//
	// c = 4 * c0;
	// c = Math.sqrt(c1 * c1 - c * c2);
	// ca = c1 + c;
	// cb = c1 - c;
	// if (Math.abs(cb) > Math.abs(ca))
	// ca = cb;
	// continue;
	// c = 2 * c0;
	// cb = c / ca;
	// x3 = x2 - cb;
	//
	// x0 = x1;
	// f0 = f1;
	// x1 = x2;
	// f1 = f2;
	// x2 = x3;
	// f2 = feval(x3); // the function using x3
	// y = Math.abs(f2);
	//
	// k = k + 1;
	//
	// } while (y > eps && k <= m);
	//
	// x = x2;
	// k = k - 1;
	//
	// }// end muller

	/**
	 * Applies the roots to the matrix
	 * 
	 * @param a
	 * @return
	 */
	public double[][] polynomial(double[][] a) {

		for (int i = 0; i < a[0].length; i++) {

			for (int j = 0; j < a.length; j++) {

				System.out.println(j);

				if (j == 0) {

					a[i][j] = Math.pow(a[i][j], 3);

				} else if (j == 1) {

					a[i][j] = Math.pow(a[i][j], 2);

				} else if (j == 2) {

					a[i][j] = Math.pow(a[i][j], 1);

				} else if (j == 3) {

					a[i][j] = Math.pow(a[i][j], 0);

				}

			}// end j for

		}// end i for
		return a;

	}// end polynomial

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

		constantLength = m.length;

		if (n[0].length > 1) {

			inverConLength = 2 * constantLength;

		} else {

			inverConLength = 1 + constantLength;

		}// end else

		// create the augmented matrix for the coefficient matrix and the
		// solution matrix

		augmentedMatrix = new double[(int) constantLength][(int) inverConLength];

		// create the augmented matrix
		// get the first matrix
		for (int i = 0; i < m.length; i++) {

			for (int j = 0; j < m[0].length; j++) {

				// add in the coefficient matrix to the augmented matrix

				augmentedMatrix[i][j] = m[i][j];

			}// end j for

		}// end i for

		double colCount = constantLength;

		// add the solutions to the the augmented matrix

		for (int i = 0; i < m.length; i++) {

			for (int j = 0; j < n[0].length; j++) {

				// add the constants to the end of the matrix
				augmentedMatrix[i][(int) (j + m[0].length)] = n[i][j];

			}// end j for

		}// end i for

		System.out.println();
		System.out.println(augmentedMatrix.length + " by "
				+ augmentedMatrix[0].length + " Matrix");

		return augmentedMatrix;

	}// end createAugmentedMatrix

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
		double determinantNum = 0;

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
	 * Make the identity matrix
	 * 
	 * @param a
	 *            the covariance matrix
	 * @return the identity of the covariance matrix
	 */
	public double[][] makeIden(double[][] a) {

		double[][] idenTemp;
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
	 * patitions the inverse from the inputed matrix
	 * 
	 * @param a
	 *            the inputed matrix
	 */
	public double[][] partitionInverse(double[][] a) {

		// initialize the matrices with the new length
		double[][] A = new double[(int) (a[0].length - constantLength)][(int) constantLength];
		double[][] inverse = new double[(int) (a[0].length - constantLength)][(int) constantLength];

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
	 * Multiplies two matrices
	 * 
	 * @param a
	 *            the first matrix
	 * @param b
	 *            the second matrix
	 * @return the product of the matrices
	 */
	public double[][] multiply(double[][] a, double[][] b) {

		System.out.println("Multiplying the two matrices");
		debug(a);
		debug(b);

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

}// end class
