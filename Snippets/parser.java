import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 * Program for lexical parser
 * 
 * @author Kavy Rattana
 *
 */
public class parser {

	boolean start;
	boolean identifier;
	boolean number;
	boolean operator;
	boolean assignment;
	boolean error;
	boolean pass;
	boolean hasSpace;
	String regex;
	String currentState;
	int[] output;

	public static void main(String args[]) {

		parser p = new parser();
		p.run();

	}// end main

	/**
	 * Method to run the methods
	 */
	public void run() {

		read();

	}// end run

	/*
	 * Reading in the text file
	 */
	public void read() {

		/*
		 * ROW 0 = start 1 = identifier 2 = number 3 = operator 4 = asssignment
		 * 5 = error
		 */

		/*
		 * COLUMN 0 = letter (a..z) 1 = digit (0..9) 2 = operator (+,-,*,/) 3 =
		 * (=) 4 = space 5 = other
		 */

		int[][] stateMachine = new int[6][6];

		// start
		stateMachine[0][0] = 1;
		stateMachine[0][1] = 2;
		stateMachine[0][2] = 3;
		stateMachine[0][3] = 4;
		stateMachine[0][4] = 0;
		stateMachine[0][5] = 5;

		// identifier
		stateMachine[1][0] = 1;
		stateMachine[1][1] = 1;
		stateMachine[1][2] = 3;
		stateMachine[1][3] = 4;
		stateMachine[1][4] = 0;
		stateMachine[1][5] = 5;

		// number
		stateMachine[2][0] = 5;
		stateMachine[2][1] = 2;
		stateMachine[2][2] = 3;
		stateMachine[2][3] = 4;
		stateMachine[2][4] = 0;
		stateMachine[2][5] = 5;

		// operator
		stateMachine[3][0] = 1;
		stateMachine[3][1] = 2;
		stateMachine[3][2] = 5;
		stateMachine[3][3] = 4;
		stateMachine[3][4] = 0;
		stateMachine[3][5] = 5;

		// assignment
		stateMachine[4][0] = 1;
		stateMachine[4][1] = 2;
		stateMachine[4][2] = 3;
		stateMachine[4][3] = 5;
		stateMachine[4][4] = 0;
		stateMachine[4][5] = 5;

		// error
		stateMachine[5][0] = 5;
		stateMachine[5][1] = 5;
		stateMachine[5][2] = 5;
		stateMachine[5][3] = 5;
		stateMachine[5][4] = 5;
		stateMachine[5][5] = 5;

		for (int i = 0; i < stateMachine.length; i++) {

			for (int j = 0; j < stateMachine[0].length; j++) {

				System.out.print("[" + stateMachine[i][j] + "]");

			}// end j for

			System.out.println();

		}// end i for

		BufferedReader br = null;

		try {

			String sCurrentLine;

			br = new BufferedReader(new FileReader("examples.txt"));

			while ((sCurrentLine = br.readLine()) != null) {

				// Parsing code start

				String[] split = sCurrentLine.split("");

				start(split, sCurrentLine, stateMachine);

			}// end while

		} catch (IOException e) {

			e.printStackTrace();

		}// end catch

	}// end read

	/**
	 * Parsing method
	 * 
	 * @param current
	 *            the current token
	 * @param s
	 *            The string
	 */
	public void start(String[] current, String s, int[][] stateMachine) {

		output = new int[current.length];
		int state = stateMachine[0][0]; // the start of the statemachine
		int stringMax = current.length; // maximum lenght of the the statement
		int i = 0, j = 0; // start state / current state
		boolean start = true;

		// establish first state
		if (current[i].matches("[a-z]")) {

			state = 1;
			output[i] = 1;

		}// end if

		if (current[i].matches("[0-9]")) {

			state = 2;
			output[i] = 2;

		}// end if

		if (current[i].equals("+") || current[i].equals("-")
				|| current[i].equals("*") || current[i].equals("/")
				&& error == false) {

			state = 3;
			output[i] = 3;

		}// end if

		if (current[i].equals("=")) {

			state = 4;
			output[i] = 4;

		}// end if

		if (current[i].equals(" ")) {

			state = 0;

		}// end if

		i++;

		// start state
		for (int k = state; k < stateMachine.length; k++) {

			for (j = 0; j < stateMachine[0].length; j++) {

				// establish first state
				if (current[i].matches("[a-z]")) {

					state = stateMachine[state][0];
					output[i] = 4;

				} else

				if (current[i].matches("[0-9]")) {

					state = stateMachine[state][1];
					output[i] = 4;

				} else

				if (current[i].equals("+") || current[i].equals("-")
						|| current[i].equals("*") || current[i].equals("/")
						&& error == false) {

					state = stateMachine[state][2];
					output[i] = 4;

				} else

				if (current[i].equals("=")) {

					state = stateMachine[state][3];
					output[i] = 4;

				} else

				if (current[i].equals(" ")) {

					state = stateMachine[state][4];
					output[i] = 4;

				} else {

					state = stateMachine[state][5];
					output[i] = 5;

				}// end else

			}// end j for

		}// end j for

		for (int m = 0; m < output.length; m++) {

			System.out.print(output[m] + " ");

		}// end m

		System.out.println();

	}// end start
}// end class
