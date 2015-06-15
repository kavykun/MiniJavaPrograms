public class ErrorHandler {

	public void nerror(int k, String string, int i, float d, String string2) {

		int c;

		System.out.println("\n*** NLIB error %i in function %s.\n    " + k
				+ " , " + string);
		switch (k) {

		case 100:
			System.out
					.println("Argument %i out of range. The value was %g. ***\n"
							+ i + ", " + d);
			break;

		case 101:
			System.out
					.println("Insufficient memory. Attempt to allocate %.0f bytes. ***\n"
							+ d);
			break;

		case 102:
			System.out.println("Could not output to file %s. ***\n" + string2);
			break;

		case 103:
			System.out.println("File not opened for input. ***\n");
			break;

		case 104:
			System.out.println("Complex division by zero. ***\n");
			break;

		case 300:
			System.out
					.println("Number of samples must be a power of two. ***\n");
			break;

		case 301:
			System.out.println("The dimension of argument %i must be odd. ***"
					+ i);
			break;

		case 302:
			System.out
					.println("The dimension of argument %i is too small. ***\n"
							+ i);
			break;

		case 303:
			System.out
					.println("Argument %i has too little frequency content. ***\n"
							+ i);
			break;

		case 304:
			System.out.println("The polynomial is not of degree %i. ***\n" + i);
			break;

		case 305:
			System.out
					.println("Arguments one and two in do not bracket a root. ***\n");
			break;

		case 306:
			System.out.println("Arguments one is not symmetric. ***\n");
			break;

		case 307:
			System.out
					.println("Maximum number of function evaluations performed. ***\n");
			break;

		case 308:
			System.out
					.println("The function F does not satisfy: F'(0) < 0. ***\n");
			break;

		case 309:
			System.out
					.println("The independent variables must be strictly increasing. ***\n");
			break;

		case 310:
			System.out.println("The initial guesses must be distinct. ***\n");
			break;

		case 311:
			System.out
					.println("Minimum step size of h = %g reached. Try increase eps\n"
							+ d);
			System.out
					.println("    or decrease length of solution interval. ***\n");
			break;

		case 312:
			System.out.println("Coefficient matrix is singular. ***\n");
			break;

		default:
			System.out.println("\nUndocumented NLIB error?\n");
			break;
		}

		/* Exit or continue ? */

		System.out.println("\nEnter x to exit, or c to continue ... ");

	}// end n error

}// end class
