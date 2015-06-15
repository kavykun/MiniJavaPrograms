import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Scanner;

/*8
 * 
 * Class to compute 1D and 2D fourier transform
 */
public class Test3part2 {

	complex complex = new complex(0, 0);
	complex[] x, result, xInv, sin, sinFFT, sinFFTinv, phase, phaseFFT,
			sinusoidalFFT, sinusoidalFFTinv, phaseFFTinv, xConjugate, limits,
			sintemp, phaseTemp, sinusoidal, sinusoidalTemp, lowFFT, lowFFTInv,
			highFFT, highFFTInv, bandFFT, bandFFTInv, notchFFT, notchFFTInv;

	int[] lowlist, highlist, band, notch;

	ErrorHandler error;
	PrintWriter pw;
	// initializing the variables
	int[] terms = { 5, 10, 1000, 2000, 5000 };
	double sampleMax = 1024, pi = Math.PI, k = 1, sum = 0.0, dx = (1 / 1024.0),
			min, max;
	String sumString;
	File file;
	Scanner s;
	int choice, currentTerm, N;

	public static void main(String[] args) {

		Test3part2 obj = new Test3part2();
		obj.run();

	}// end main

	/**
	 * Method to run the program
	 */
	public void run() {

		do {

			System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
			System.out.println("Choose a operation to run");
			System.out.println("Run FFT-1 with the coresponding function");
			System.out.println("1. Run 2a.i fs on S = 5");
			System.out.println("2. Run 2a.i fs on S = 10");
			System.out.println("3. Run 2a.i fs on S = 1000");
			System.out.println("4. Run 2a.i gs on S = 5");
			System.out.println("5. Run 2a.i gs on S = 10");
			System.out.println("6. Run 2a.i gs on S = 1000");
			System.out.println("7. 2a.ii FFT f");
			System.out.println("8. 2a.ii FFT g");
			System.out.println("9. 2a.ii FFT f inverse");
			System.out.println("10. 2a.ii FFT g inverse");
			System.out.println("11. 2a.ii FFTf*FFTf-1");
			System.out.println("12. 2a.ii FFTg*FFTg-1");
			System.out.println("13. Print limitsf");
			System.out.println("14. Print limitsg");
			System.out.println("15. Sin Sum");
			System.out.println("16. Sin Product");
			System.out.println("17. Phase");
			System.out.println("18. Sinusoidal tone");
			System.out.println("19. Low-pass filter");
			System.out.println("20. High=pass filter");
			System.out.println("21. Band-pass");
			System.out.println("22. Notch filter");
			System.out.println("23. Correlation");
			System.out.println("24. 2D FFT");
			System.out.println("25. PSD of limit functions");
			System.out.println("26. Exit Program");
			System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");

			s = new Scanner(System.in);
			choice = s.nextInt();

			switch (choice) {

			// S = 5
			case 1:
				file = new File("Sf5.txt");
				N = 1024;
				x = new complex[N];

				try {

					pw = new PrintWriter(new FileWriter(file), true);
					currentTerm = 0;
					x = runAlg2f(x);
					printToFile(x);

				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				pw.close();
				break;

			// S = 10
			case 2:
				file = new File("Sf10.txt");
				N = 1024;
				x = new complex[N];

				try {

					pw = new PrintWriter(new FileWriter(file), true);
					currentTerm = 1;
					x = runAlg2f(x);
					printToFile(x);

				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				pw.close();
				break;

			// S = 1000
			case 3:
				file = new File("Sf1000.txt");
				N = 1024;
				x = new complex[N];

				try {

					pw = new PrintWriter(new FileWriter(file), true);
					currentTerm = 2;
					x = runAlg2f(x);
					printToFile(x);

				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				pw.close();
				break;

			// S = 5
			case 4:
				file = new File("Sg5.txt");
				N = 1024;
				x = new complex[N];

				try {

					pw = new PrintWriter(new FileWriter(file), true);
					currentTerm = 0;
					x = runAlg2g(x);
					printToFile(x);

				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				pw.close();
				break;

			// S = 10
			case 5:
				file = new File("Sg10.txt");
				N = 1024;
				x = new complex[N];

				try {

					pw = new PrintWriter(new FileWriter(file), true);
					currentTerm = 1;
					x = runAlg2g(x);
					printToFile(x);

				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				pw.close();
				break;

			// S = 1000
			case 6:
				file = new File("Sg1000.txt");
				N = 1024;
				x = new complex[N];

				try {

					pw = new PrintWriter(new FileWriter(file), true);
					currentTerm = 4;
					x = runAlg2g(x);
					printToFile(x);

				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				pw.close();
				break;

			// FFT of ff
			case 7:
				file = new File("complexf.txt");

				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				x = new complex[N];

				// original data
				currentTerm = 4;
				x = runAlg2f(x);

				// x[0] = new complex((float) 26160.0, 0);
				// x[1] = new complex((float) 19011.0, 0);
				// x[2] = new complex((float) 18757.0, 0);
				// x[3] = new complex((float) 18405.0, 0);
				// x[4] = new complex((float) 17888.0, 0);
				// x[5] = new complex((float) 14720.0, 0);
				// x[6] = new complex((float) 14285.0, 0);
				// x[7] = new complex((float) 17018.0, 0);
				// x[8] = new complex((float) 18014.0, 0);
				// x[9] = new complex((float) 17119.0, 0);
				// x[10] = new complex((float) 16400.0, 0);
				// x[11] = new complex((float) 17497.0, 0);
				// x[12] = new complex((float) 17846.0, 0);
				// x[13] = new complex((float) 15700.0, 0);
				// x[14] = new complex((float) 17636.0, 0);
				// x[15] = new complex((float) 17181.0, 0);

				System.out.println("\nThe Data");
				for (int i = 0; i < x.length; i++) {

					System.out.println(x[i].x);

				}// end i for

				System.out.println("\nThe Complex Matrix");
				printComplexArray(x);

				// run the FFT

				result = FFT(x, N, 1);
				printToFile(result);

				pw.close();
				break;

			// FFT of fg
			case 8:
				file = new File("complexg.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				x = new complex[N];

				// original data
				currentTerm = 4;
				x = runAlg2g(x);

				System.out.println("\nThe Data");
				for (int i = 0; i < x.length; i++) {

					System.out.println(x[i].x);

				}// end i for

				System.out.println("\nThe Complex Matrix");
				printComplexArray(x);

				// run the FFT
				result = FFT(x, N, 1);
				printToFile(result);

				pw.close();

				break;

			// inverse of fftf
			case 9:

				file = new File("complexf-1.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				xInv = new complex[N];
				xConjugate = new complex[N];

				// original data
				currentTerm = 4;
				xInv = runAlg2f(xInv);

				System.out.println("\nThe Data");
				for (int i = 0; i < xInv.length; i++) {

					System.out.println(xInv[i].x);

				}// end i for

				// run the FFT
				result = FFT(xInv, N, 1);

				xConjugate = conjugate(result);

				printToFile(xConjugate);

				printComplexArray(xConjugate);

				pw.close();

				break;

			// inverse of fftg
			case 10:

				file = new File("complexg-1.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				xInv = new complex[N];
				xConjugate = new complex[N];

				// original data
				currentTerm = 4;
				xInv = runAlg2g(xInv);

				System.out.println("\nThe Data");
				for (int i = 0; i < xInv.length; i++) {

					System.out.println(xInv[i].x);

				}// end i for

				// run the FFT
				result = FFT(xInv, N, 1);

				xConjugate = conjugate(result);

				printToFile(xConjugate);

				printComplexArray(xConjugate);

				pw.close();

				break;
			// product of fftf and fftfinv
			case 11:

				file = new File("complexffinv.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				result = multiplyFFT(x, xConjugate);
				printToFile(result);

				pw.close();
				break;

			// product of fftg and fftginv
			case 12:

				file = new File("complexgginv.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				result = multiplyFFT(x, xConjugate);
				printToFile(result);

				pw.close();

				break;

			// limits functions f
			case 13:

				file = new File("limitsfmin.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				result = new complex[N];

				currentTerm = 2;

				limits = limitsf(result, min);

				printToFile(limits);

				file = new File("limitsfmax.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				N = 1024;
				result = new complex[N];

				currentTerm = 2;

				limits = limitsf(result, max);

				printToFile(limits);

				pw.close();

				break;

			// limits functions g
			case 14:

				file = new File("limitsgmin.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				result = new complex[N];

				currentTerm = 2;

				limits = limitsg(result, min);

				printToFile(limits);

				file = new File("limitsgmax.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				N = 1024;
				result = new complex[N];

				currentTerm = 2;

				limits = limitsg(result, max);

				printToFile(limits);

				pw.close();

				break;
			// sin sum
			case 15:

				file = new File("sinsum.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				N = 1024;
				sin = new complex[N];

				sintemp = sin(sin);

				sinFFT = FFT(sintemp, N, 1);
				sinFFTinv = conjugate(sinFFT);

				result = multiplyFFT(sinFFT, sinFFTinv);
				printToFile(result);

				pw.close();

				break;

			// sin product
			case 16:

				file = new File("sinproduct.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				N = 1024;
				sin = new complex[N];

				sintemp = sin(sin);

				sinFFT = FFT(sintemp, N, 1);
				sinFFTinv = conjugate(sinFFT);

				result = multiplyFFT(sinFFT, sinFFTinv);
				printToFile(result);

				pw.close();

				break;

			// phase
			case 17:

				file = new File("phase.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				N = 1024;
				phase = new complex[N];

				phaseTemp = phase(phase);

				// printComplexArray(phaseTemp);

				phaseFFT = FFT(phaseTemp, N, 1);
				phaseFFTinv = conjugate(phaseFFT);

				result = multiplyFFT(phaseFFT, phaseFFTinv);
				printToFile(result);

				pw.close();

				break;

			// sinusoidal
			case 18:

				file = new File("sinusoidal.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				N = 1024;
				sinusoidal = new complex[N];

				sinusoidalTemp = sinusoidal(sinusoidal);

				sinusoidalFFT = FFT(sinusoidalTemp, N, 1);

				System.out.println("working");
				printToFile(sinusoidalFFT);

				// sinusoidalFFTinv = conjugate(sinusoidalFFT);
				//
				// result = multiplyFFT(sinusoidalFFT, sinusoidalFFTinv);
				// printToFile(result);

				pw.close();

				break;

			// low pass
			case 19:

				file = new File("lowpass.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				x = new complex[N];
				lowlist = new int[5];

				// original data
				currentTerm = 2;
				x = runAlg2f(x);

				System.out.println("\nThe Data");
				for (int i = 0; i < x.length; i++) {

					System.out.println(x[i].x);

				}// end i for

				lowFFT = FFT(x, N, 1);

				System.out.println("\nThe Complex Matrix");
				printComplexArray(x);

				complex[] lowFFTTemp = new complex[lowFFT.length];

				lowFFTTemp = lowpass(lowFFT); // pass through filter

				lowFFTInv = FFT(lowFFTTemp, N, -1);

				printToFile(lowFFTInv);

				break;

			// high pass
			case 20:

				file = new File("highpass.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				x = new complex[N];

				// original data
				currentTerm = 2;
				x = runAlg2f(x);

				System.out.println("\nThe Data");
				for (int i = 0; i < x.length; i++) {

					System.out.println(x[i].x);

				}// end i for

				highFFT = FFT(x, N, 1);

				System.out.println("\nThe Complex Matrix");
				printComplexArray(x);

				complex[] highFFTTemp = new complex[highFFT.length];

				highFFTTemp = highpass(highFFT); // pass through filter

				highFFTInv = FFT(highFFTTemp, N, -1);

				printToFile(highFFTInv);

				break;

			case 21:

				file = new File("bandpass.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				x = new complex[N];

				// original data
				currentTerm = 2;
				x = runAlg2f(x);

				System.out.println("\nThe Data");
				for (int i = 0; i < x.length; i++) {

					System.out.println(x[i].x);

				}// end i for

				bandFFT = FFT(x, N, 1);

				System.out.println("\nThe Complex Matrix");
				printComplexArray(x);

				complex[] bandFFTTemp = new complex[bandFFT.length];

				bandFFTTemp = bandpass(bandFFT); // pass through filter

				bandFFTInv = FFT(bandFFTTemp, N, -1);

				printToFile(bandFFTInv);

				break;

			// notch
			case 22:

				file = new File("notchpass.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				N = 1024;
				x = new complex[N];

				// original data
				currentTerm = 2;
				x = runAlg2f(x);

				System.out.println("\nThe Data");
				for (int i = 0; i < x.length; i++) {

					System.out.println(x[i].x);

				}// end i for

				notchFFT = FFT(x, N, 1);

				System.out.println("\nThe Complex Matrix");
				printComplexArray(x);

				complex[] notchFFTTemp = new complex[notchFFT.length];

				notchFFTTemp = notchpass(notchFFT); // pass through filter

				notchFFTInv = FFT(notchFFTTemp, N, -1);

				printToFile(notchFFTInv);

				break;

			// correlation
			case 23:

				// Assume that the receiver is turned off for 0.04 seconds after
				// the pulse is transmitted, 100 kHz sampling rate, and the
				// return signal measurements were taken from the first
				// 1024-sample window after the receiver begins listening again.
				// How far away is the primary reflector (6 points)?

				String input;

				file = new File("correlation.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				N = 1024;

				complex[] pulse = new complex[N];
				complex[] signal = new complex[N];

				complex[] pulseFFT = new complex[N];
				complex[] signalFFT = new complex[N];

				complex[] pulseFFTInv = new complex[N];
				complex[] output = new complex[N];
				complex[] product = new complex[N];

				complex[] pulseData = new complex[N];
				complex[] signalData = new complex[N];

				complex[] plot = new complex[N / 4];

				pulse = new complex[N];
				signal = new complex[N];

				Scanner s = new Scanner(System.in);
				System.out.println("Enter the file name for the pulse: ");
				input = s.next();
				pulse = readFile(pulseData, input);

				printComplexArray(pulse);

				System.out.println("Enter the file name for the signal: ");
				input = s.next();
				signal = readFile(signalData, input);

				printComplexArray(signal);

				file = new File("signalPlot.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				printToFile(signal);

				// FFT of the signal
				// FFT of the pulse
				// zeros into the rest of the pulse
				// FFT of signal * conjugate of the pulse
				// FFT inverse of that

				signalFFT = FFT(signal, N, 1);
				pulseFFT = FFT(pulse, N, 1);
				pulseFFTInv = conjugate(pulseFFT);

				product = multiplyFFT(signalFFT, pulseFFTInv);

				output = FFT(product, N, -1);

				printToFile(output);

				// convolation

				// take one quater of the results - fist 256 values
				// make anotehr complex arrray
				// p = 20
				// make array of complex values where
				// complex number 1/a a = passsed , 0
				// first 19 values , rest if zero

				// FFT of both of them
				// Multiply them together
				// FFT inverse

				file = new File("convoltion.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				N = 1024;

				complex[] convolution = new complex[256];
				complex[] pfilter = new complex[256];

				complex[] convolutionFFT = new complex[256];
				complex[] pfilterFFT = new complex[256];

				complex[] convolutionInv = new complex[256];

				float a = 20;

				System.out.println("Enter the file name for the signal: ");
				input = s.next();
				signal = readFile(signalData, input);

				printComplexArray(signal);

				pw.close();

				file = new File("signalPlotOverlay.txt");
				try {
					pw = new PrintWriter(new FileWriter(file), true);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				int index = 0;

				for (int i = 0; i < signal.length; i++) {

					if (i >= 97 && i < 353) {

						convolution[index] = signal[i];
						plot[index] = signal[i];
						index++;

					}// end if

				}// end i for

				printToFile(plot);

				for (int j = 0; j < pfilter.length; j++) {

					if (j < 20) {

						pfilter[j] = new complex((float) 1 / a, 0);

					} else {

						pfilter[j] = new complex(0, 0);

					}// end else

				}// end j for

				pfilterFFT = FFT(pfilter, 256, 1);
				convolutionFFT = FFT(convolution, 256, 1);

				product = multiplyFFT(convolutionFFT, pfilterFFT);

				convolutionInv = FFT(product, 256, -1);

				printToFile(convolutionInv);

				pw.close();

				break;

			case 24:

				complex[][] twod = new complex[512][512];

				TwoDFFT();

				break;

			case 25:

				s = new Scanner(System.in);

				complex[] f5,
				f10,
				f1000,
				g5,
				g10,
				g1000,
				psd;

				// F

				// System.out.println("Enter the file name for the f5: ");
				// input = s.next();
				// signal = readFile(signalData, input);
				//
				// System.out.println("Enter the file name for the f10: ");
				//
				// input = s.next();
				// signal = readFile(signalData, input);
				//
				// System.out.println("Enter the file name for the f1000: ");
				//
				// input = s.next();
				// signal = readFile(signalData, input);
				//
				// // G
				//
				// System.out.println("Enter the file name for the f5: ");
				// input = s.next();
				// signal = readFile(signalData, input);
				//
				// System.out.println("Enter the file name for the f10: ");
				//
				// input = s.next();
				// signal = readFile(signalData, input);
				//
				// System.out.println("Enter the file name for the f1000: ");
				//
				// input = s.next();
				// signal = readFile(signalData, input);

				// printComplexArray(signal);

				break;

			default:
				System.err.println("Invalid Choice. Try again");
				break;

			}// end switch

		} while (choice >= 1 && choice <= 30);

	}// end run

	/**
	 * Runs the algorithm for 2a.iii Running the algorithm on sample size =
	 * 1000;
	 */
	public complex[] runAlg2f(complex[] a) {

		max = -99999999;
		min = 99999999;

		// running the algorithms
		for (int i = 0; i < sampleMax; i++) {

			double tempi = i * dx;

			for (k = 1; k <= terms[currentTerm]; k++) {

				sum += (Math.sin((2 * pi) * (2 * k - 1) * tempi) / (2 * k - 1));

			}// end j for

			if (sum > max) {

				max = sum;

			}

			if (sum < min) {

				min = sum;

			}

			a[i] = new complex((float) sum, 0);

			sum = 0.0;

		}// end i for

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		System.out.println("length: " + a.length);
		System.out.println("first number: " + a[0].x);
		System.out.println("last number: " + a[a.length - 1].x);
		System.out.println("Max: " + max);
		System.out.println("Max: " + min);

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		return a;

	}// end runAlg3

	/**
	 * Runs the algorithm for 2a.iii Running the algorithm on size sample size =
	 * 1000;
	 */
	public complex[] runAlg2g(complex[] a) {

		max = -99999999;
		min = 99999999;

		// running the algorithms
		for (int i = 0; i < sampleMax; i++) {

			double tempi = i * dx;

			for (k = 1; k <= terms[currentTerm]; k++) {

				sum += (Math.sin((2 * pi) * (2 * k) * tempi) / (2 * k));

			}// end j for

			if (sum > max) {

				max = sum;

			}

			a[i] = new complex((float) sum, 0);

			sum = 0.0;

		}// end i for

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		System.out.println("length: " + a.length);
		System.out.println("first number: " + a[0].x);
		System.out.println("last number: " + a[a.length - 1].x);
		System.out.println("Max: " + max);
		System.out.println("Min " + min);

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		return a;

	}// end runAlg3

	/**
	 * Prints the limits
	 * 
	 */
	public complex[] limitsf(complex[] a, double n) {

		// running the algorithms
		for (int i = 0; i < a.length; i++) {

			a[i] = new complex((float) n, 0);

		}// end i for

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		System.out.println("length: " + a.length);
		System.out.println("first number: " + a[0].x);
		System.out.println("last number: " + a[a.length - 1].x);
		System.out.println("Max: " + max);
		System.out.println("Max: " + min);

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		return a;

	}// end limitsf

	/**
	 * Prints the limits
	 * 
	 */
	public complex[] limitsg(complex[] a, double n) {

		// running the algorithms
		for (int i = 0; i < a.length; i++) {

			a[i] = new complex((float) n, 0);

		}// end i for

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		System.out.println("length: " + a.length);
		System.out.println("first number: " + a[0].x);
		System.out.println("last number: " + a[a.length - 1].x);
		System.out.println("Max: " + max);
		System.out.println("Max: " + min);

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		return a;

	}// end limits

	/**
	 * Prints the complex[] to a file
	 * 
	 * @param a
	 */
	public void printToFile(complex[] a) {

		for (int i = 0; i < a.length; i++) {

			if (choice == 11 || choice == 12 || choice == 15 || choice == 16) {

				if (a[i].x >= 0 && i < (a.length / 2)) {

					pw.println(a[i].x + "");

				}// end if

			} else {

				pw.println(a[i].x);

			}// end if

		}// end i for

	}// end printToFile

	/**
	 * Prints the complex[] to a file
	 * 
	 * @param a
	 */
	public void printToFiletwod(complex[][] a) {

		complex[][] z = new complex[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				z[i][j] = a[i][j];

			}// end j for

		}// end i for

		for (int i = 0; i < z.length; i++) {

			for (int j = 0; j < z[0].length; j++) {

				pw.print("[" + z[i][j].x + " " + z[i][j].y + "]");

			}// end if

			pw.println();

		}// end i for

	}// end printToFile

	/**
	 * Method to run the FFT
	 */
	public complex[] FFT(complex[] a, int n, int d) {

		// complex[] z = new complex[a.length];

		complex[] z = new complex[a.length];

		for (int i = 0; i < a.length; i++) {

			z[i] = a[i];

		}// end i for

		int i = 1, k, m, r, i0, i1;
		float PI = (float) Math.PI;
		complex w, u, t;
		float theta;

		// for (int h = 0; h < a.length; h++) {
		//
		// z[h] = a[h];
		//
		// }

		/* Initialize */

		if (n < 2)
			error.nerror(100, "dft", 2, (float) n, "");
		else if ((d != -1) && (d != 1))
			error.nerror(100, "dft", 3, (float) d, "");
		while (i < n)
			i *= 2;
		if (i > n)
			error.nerror(300, "dft", 0, 0, "");
		theta = -d * 2 * PI / n;
		r = n / 2;

		for (i = 1; i < n - 1; i *= 2) {
			w = complex.cnum((float) Math.cos(i * theta),
					(float) Math.sin(i * theta));
			for (k = 0; k < n - 1; k += 2 * r) {
				u = complex.cnum(1, 0);
				for (m = 0; m < r; m++) {
					i0 = k + m;
					i1 = i0 + r;
					t = complex.csub(z[i0], z[i1]);
					z[i0] = complex.cadd(z[i0], z[i1]);
					z[i1] = complex.cmul(t, u);
					u = complex.cmul(u, w);

					// System.out.println("i=" + i + "; " + "k=" + k + "; m=" +
					// m
					// + "; r=" + r);
				}
			}
			r /= 2;
		}

		// System.out.println("Results from first loop: ");
		// printComplexArray(z);

		// System.out.println("Rearranging the results");
		for (i = 0; i < n - 1; i++) {
			r = i;
			k = 0;
			for (m = 1; m < n - 1; m *= 2) {
				k = 2 * k + (r % 2);
				r /= 2;

				// System.out.println("i=" + i + "; " + "k=" + k + "; m=" + m
				// + "; r=" + r);
			}
			if (k > i) {
				t = z[i];
				z[i] = z[k];
				z[k] = t;
			}
		}

		if (d < 0) {
			for (i = 1; i <= n - 1; i++) {
				z[i].x /= n;
				z[i].y /= n;
			}
		}

		// System.out.println("The FFT of the Data");
		// printComplexArray(z);

		return z;

	}// end FFT

	/**
	 * Multiplies two FFTs together
	 * 
	 * @return the product of the two complex array
	 */
	public complex[] multiplyFFT(complex[] a, complex[] b) {

		complex[] output = new complex[a.length];

		if (a.length != b.length) {

			System.err.println("The two complex arrays do not match");

		} else {

			for (int i = 0; i < a.length; i++) {

				output[i] = complex.cmul(a[i], b[i]);

			}// end i for

		}// end else

		return output;

	}

	/**
	 * Multiplies two 2DFFTs together
	 * 
	 * @return the product of the two complex array
	 */
	public complex[][] multiplyFFTTwoD(complex[][] a, complex[][] b) {

		complex[][] output = new complex[a.length][a[0].length];

		complex[][] first = new complex[a.length][a[0].length];
		complex[][] second = new complex[b.length][b[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a.length; j++) {

				first[i][j] = a[i][j];

			}// end j for
		}// end i for

		for (int i = 0; i < b.length; i++) {

			for (int j = 0; j < b.length; j++) {

				second[i][j] = b[i][j];

			}// end j for
		}// end i fo
		if (a.length != b.length) {

			System.err.println("The two complex arrays do not match");

		} else {

			for (int i = 0; i < a.length; i++) {

				for (int j = 0; j < a[0].length; j++) {

					second[i][j].negateY();
					output[i][j] = complex.cmul(first[i][j], second[i][j]);

				}// end j for

			}// end i for

		}// end else

		return output;

	}// end multiplyFFTTwoD

	/**
	 * Adds two FFTs together
	 * 
	 * @return the resulting sum of complex array
	 */
	public complex[] addFFT(complex[] a, complex[] b) {

		complex[] output = new complex[a.length];

		if (a.length != b.length) {

			System.err.println("The two complex arrays do not match");

		} else {

			for (int i = 0; i < a.length; i++) {

				output[i] = complex.cadd(a[i], b[i]);

			}// end i for

			System.out.println("Printing complex: ");
			printComplexArray(output);

		}// end else

		return output;

	}

	/**
	 * Runs the sin functions Sum and products
	 * 
	 * @param z
	 *            the complex array to be evaluated
	 * @return a the resulting complex array
	 */
	public complex[] sin(complex[] z) {

		complex[] x = new complex[z.length];
		complex[] y = new complex[z.length];
		double product1 = 0, product2 = 0, product = 0;
		double a, f1, f2, t, c, PI;
		PI = Math.PI;
		a = 1;
		c = 0;
		f1 = 10;
		f2 = 43;

		// running the algorithms
		for (int i = 0; i < sampleMax; i++) {

			double tempi = i * dx;

			if (choice == 15) {

				System.out.println(i);
				product1 += a * Math.sin(2 * PI * f1 * (tempi - c));
				product2 += a * Math.sin(2 * PI * f2 * (tempi - c));

			} else if (choice == 16) {

				product1 += a * Math.sin(2 * PI * f1 * (tempi - c));
				product2 += a * Math.sin(2 * PI * f2 * (tempi - c));

			}// end else if13

			x[i] = new complex((float) product1, 0);
			y[i] = new complex((float) product2, 0);

			product1 = 0.0;
			product2 = 0.0;
			product = 0.0;

		}// end i for

		if (choice == 15) {

			z = addFFT(y, x);

		} else if (choice == 16) {

			z = multiplyFFT(x, y);

		}// end else if

		return z;

	}// end sin

	/**
	 * The phase method
	 * 
	 * @param z
	 *            the complex array to be evaluated
	 * @return the resulting complex array
	 */
	public complex[] phase(complex[] z) {

		double max = -99999999;

		int rand = 0 + (int) (Math.random() * ((1024 - 0) + 1));

		System.out.println(rand);

		// running the algorithms
		for (int i = 0; i < sampleMax; i++) {

			if (i == rand) {

				z[i] = new complex((float) 1, 0);

			} else {

				z[i] = new complex((float) 0, 0);

			}// end else

		}// end i for

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		System.out.println("length: " + z.length);
		System.out.println("first number: " + z[0].x);
		System.out.println("last number: " + z[z.length - 1].x);
		System.out.println("Max: " + max);

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		return z;

	}// end phase

	public complex[] sinusoidal(complex[] z) {

		double max = -99999999;
		double PI = Math.PI, product;
		double c = 0.1;

		// running the algorithms
		for (int i = 0; i < sampleMax; i++) {

			double tempi = i * dx;

			product = Math.sin(38 * PI * (tempi - c));

			z[i] = new complex((float) product, 0);

		}// end i for

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		System.out.println("length: " + z.length);
		System.out.println("first number: " + z[0].x);
		System.out.println("last number: " + z[z.length - 1].x);
		System.out.println("Max: " + max);

		System.out.println("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		return z;

	}

	/**
	 * Runs a low pass filter with the signal
	 * 
	 * @param z
	 *            the signal
	 * @return the altered signal
	 */
	public complex[] lowpass(complex[] z) {

		complex[] output = new complex[z.length];
		complex one = new complex(1, 0);
		complex zero = new complex(0, 0);

		for (int i = 0; i < z.length; i++) {

			if (i == 1 || i == 3 || i == 5 || i == 7 || i == 9 || i == 1023
					|| i == 1021 || i == 1019 || i == 1017 || i == 1015) {

				output[i] = complex.cmul(z[i], one);

			} else {

				output[i] = complex.cmul(z[i], zero);

			}// end else

		}// end i for

		return output;

	}// end low pass

	/**
	 * Runs a low pass filter with the signal
	 * 
	 * @param z
	 *            the signal
	 * @return the altered signal
	 */
	public complex[] highpass(complex[] z) {

		complex[] output = new complex[z.length];
		complex one = new complex(1, 0);
		complex zero = new complex(0, 0);

		for (int i = 0; i < z.length; i++) {

			if (i == 1 || i == 3 || i == 5 || i == 7 || i == 9 || i == 1023
					|| i == 1021 || i == 1019 || i == 1017 || i == 1015) {

				output[i] = complex.cmul(z[i], zero);

			} else {

				output[i] = complex.cmul(z[i], one);

			}// end else

		}// end i for

		return output;

	}

	/**
	 * Runs a low pass filter with the signal
	 * 
	 * @param z
	 *            the signal
	 * @return the altered signal
	 */
	public complex[] bandpass(complex[] z) {

		complex[] output = new complex[z.length];
		complex one = new complex(1, 0);
		complex zero = new complex(0, 0);

		for (int i = 0; i < z.length; i++) {

			if (i == 4 || i == 5 || i == 6 || i == 7 || i == 1021 || i == 1020
					|| i == 1019 || i == 1018) {

				output[i] = complex.cmul(z[i], one);

			} else {

				output[i] = complex.cmul(z[i], zero);

			}// end else

		}// end i for

		return output;
	}

	/**
	 * Runs a low pass filter with the signal
	 * 
	 * @param z
	 *            the signal
	 * @return the altered signal
	 */
	public complex[] notchpass(complex[] z) {

		complex[] output = new complex[z.length];
		complex one = new complex(1, 0);
		complex zero = new complex(0, 0);

		for (int i = 0; i < z.length; i++) {

			if (i == 4 || i == 5 || i == 6 || i == 7 || i == 1021 || i == 1020
					|| i == 1019 || i == 1018) {

				output[i] = complex.cmul(z[i], zero);

			} else {

				output[i] = complex.cmul(z[i], one);

			}// end else

		}// end i for

		return output;
	}

	/**
	 * Runs the 2D FFT with the complex[][]
	 * 
	 * @return
	 */
	public complex[][] TwoDFFT() {

		complex[][] pixelsSignal = new complex[512][512];
		complex[][] pixelsPulse = new complex[512][512];

		N = 512;
		Color white = new Color(255, 255, 255);

		Picture picSignal = new Picture(N, N);
		Picture picOutput = new Picture(N, N);

		for (int i = 230; i < 230 + 150; i++) {

			for (int j = 260; j < 260 + 75; j++) {

				picSignal.set(i, j, white);

			}// end j for

		}// end i for

		picSignal.show();

		Picture picPulse = new Picture(N, N);

		for (int i = 0; i < 50; i++) {

			for (int j = 0; j < 25; j++) {

				picPulse.set(i, j, white);

			}// end j for

		}// end i for

		picPulse.show();

		// get the signal colors

		Color[][] colorArray = new Color[512][512];

		colorArray = picSignal.getColorArray();

		for (int m = 0; m < colorArray.length; m++) {

			for (int l = 0; l < colorArray[0].length; l++) {

				pixelsSignal[m][l] = new complex(colorArray[m][l].getRed(), 0);

			}// end l

		}// end m

		// get the pusle colors

		colorArray = picPulse.getColorArray();

		for (int m = 0; m < colorArray.length; m++) {

			for (int l = 0; l < colorArray[0].length; l++) {

				pixelsPulse[m][l] = new complex(colorArray[m][l].getRed(), 0);

			}// end l

		}// end m

		/* Compute m FFTs of length n */

		complex[][] pixelsSignalFFT = new complex[512][512];
		complex[][] pixelsPulseFFT = new complex[512][512];
		complex[][] pixelsConjugate = new complex[512][512];
		complex[][] productTwoD = new complex[512][512];
		complex[][] productInv = new complex[512][512];

		complex[][] pixelsPulseTemp = new complex[512][512];

		for (int i = 0; i < pixelsPulseTemp.length; i++) {

			for (int j = 0; j < pixelsPulseTemp[0].length; j++) {

				pixelsPulseTemp[i][j] = pixelsPulse[i][j];

			}// end j for

		}// end i for
			// //////////////////////////////////////////////////////Signal////////////////////////////

		pixelsSignalFFT = algorithm(pixelsSignal, 1);
		pixelsPulseFFT = algorithm(pixelsPulse, 1);
		pixelsConjugate = conjugateTwoD(pixelsPulseFFT);

		System.out.println("Pulse");
		// printtwod(pixelsPulse);
		System.out.println("Signal");
		// printtwod(pixelsSignal);

		//
		// The conjugate of the pulseFFT

		// Multiply the signal with the pulse conjugate
		productTwoD = multiplyFFTTwoD(pixelsSignalFFT, pixelsConjugate);

		// ////////////////////////////////////////////////////////////////////////////////////
		// Take the inverse of the product
		productInv = algorithm(productTwoD, -1);

		// Find the max of the inverse

		double max = -999999999;
		double min = 999999999;

		for (int i = 0; i < productInv.length; i++) {

			for (int j = 0; j < productInv[0].length; j++) {

				if (productInv[i][j].x > max) {

					max = productInv[i][j].x;

				}// end if

			}// end j for

		}// end i for

		// ////////////////////////////////ok/////////////////////

		for (int i = 0; i < productInv.length; i++) {

			for (int j = 0; j < productInv[0].length; j++) {

				if (productInv[i][j].x < min) {

					min = productInv[i][j].x;

				}// end if

			}// end j
		}// end i

		printtwod(productInv);

		System.out.println("max = " + max);
		System.out.println("min = " + min);

		// double mx = Math.log(max - min + 2);
		double mx = (255 / max);
		double out;

		// ////////////////////////////////////////////////////

		for (int i = 0; i < productInv.length; i++) {

			for (int j = 0; j < productInv[0].length; j++) {

				// out = Math.log(productInv[i][j].x - min + 1) * (255 / mx);
				out = productInv[i][j].x * mx;

				productInv[i][j].setX((int) out);

			}// end j for

		}// end i for

		printtwod(productInv);

		Picture resultingPicture = new Picture(512, 512);

		for (int i = 0; i < resultingPicture.height(); i++) {

			for (int j = 0; j < resultingPicture.width(); j++) {

				Color color = new Color((int) productInv[i][j].getX(),
						(int) productInv[i][j].getX(),
						(int) productInv[i][j].getX());

				resultingPicture.set(i, j, color);

			}// end j for

		}// end i for

		resultingPicture.show();

		Picture redPicture = new Picture(512, 512);

		for (int i = 0; i < redPicture.height(); i++) {

			for (int j = 0; j < redPicture.width(); j++) {

				if (productInv[i][j].getX() >= 229.5
						&& productInv[i][j].getX() < 255) {

					Color red = new Color(255, 0, 0);
					redPicture.set(i, j, red);

				} else {

					Color color = new Color((int) productInv[i][j].getX(),
							(int) productInv[i][j].getX(),
							(int) productInv[i][j].getX());

					redPicture.set(i, j, color);

				}// end else

			}// end j for
		}// end i for

		redPicture.show();

		// printtwod(product);
		//
		// for (i = 0; i < product.length; i++) {
		//
		// for (int j = 0; j < product[0].length; j++) {
		//
		// if (product[i][j].x == 255) {
		//
		// picOutput.set(i, j, Color.WHITE);
		//
		// } else {
		//
		// picOutput.set(i, j, Color.BLACK);
		//
		// }

		// }// end j for
		//
		// }// end i for

		// System.out.println("Showing output");
		// picOutput.show();

		return null;

	}// end TwoDFFT

	/**
	 * Runs the algorithm for the 2D FFT
	 * 
	 * @return
	 */
	public complex[][] algorithm(complex[][] pixelsArray, int b) {

		// //////////////////////////////////////////////////////Signal///////////////////////////

		int d = b;

		complex[][] Y;
		complex[][] Z = new complex[pixelsArray.length][pixelsArray[0].length];

		for (int i = 0; i < pixelsArray.length; i++) {

			for (int j = 0; j < pixelsArray[0].length; j++) {

				Z[i][j] = pixelsArray[i][j];

			}// end j for

		}// end i for
		int m = pixelsArray.length, n = pixelsArray[0].length;

		/* Initialize */

		complex[] row1FFT = new complex[pixelsArray.length];
		complex[] col1FFT = new complex[pixelsArray.length];
		complex[] row = new complex[pixelsArray.length];
		complex[] col = new complex[pixelsArray.length];
		Y = new complex[pixelsArray.length][pixelsArray[0].length];

		int i, j, k;
		/* Compute m FFTs of length n */

		/* Compute m FFTs of length n */

		for (k = 0; k < m; k++) {
			for (i = 0; i < n; i++)
				row[i] = Z[k][i];
			row1FFT = FFT(row, n, d);
			for (i = 0; i < n; i++)
				Y[i][k] = row1FFT[i];
		}

		/* Compute n FFTs of length m */

		for (i = 0; i < n; i++) {
			for (k = 0; k < m; k++)
				col[k] = Y[i][k];
			col1FFT = FFT(col, m, d);
			for (k = 0; k < m; k++)
				Y[i][k] = col1FFT[k];
		}

		/* Finalize */

		for (k = 0; k < m; k++)
			for (i = 0; i < n; i++)
				Z[k][i] = Y[i][k];

		return Z;

	}

	/**
	 * Prints out the complex array
	 */
	public void printComplexArray(complex[] a) {

		for (int h = 0; h < a.length; h++) {

			System.out.println(a[h].x + "  " + a[h].y + "i");

		}// end i for

	}// end printComplexArray

	/**
	 * Find the conjugate of the complex array
	 * 
	 * @param z
	 *            the complex array
	 * @return the conjugated complex array
	 */
	public complex[] conjugate(complex[] z) {

		for (int i = 0; i < z.length; i++) {

			double conjugate = z[i].y * -1;
			z[i].setY(conjugate);

		}// end i for

		return z;

	}// end conjugate

	/**
	 * Find the conjugate of the complex array
	 * 
	 * @param z
	 *            the complex array
	 * @return the conjugated complex array
	 */
	public complex[][] conjugateTwoD(complex[][] a) {

		complex[][] z = new complex[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				z[i][j] = a[i][j];

			}// end j for

		}// end i for

		for (int i = 0; i < z.length; i++) {

			for (int j = 0; j < z[0].length; j++) {

				z[i][j].negateY();

			}// end j for

		}// end i for

		return z;

	}// end conjugate

	/**
	 * Reads in the file
	 */
	public complex[] readFile(complex[] z, String f) {

		File file = new File(f);
		String line;
		FileReader fr;
		BufferedReader br = null;
		int i = 0;

		try {

			fr = new FileReader(file);
			br = new BufferedReader(fr);

			try {

				while ((line = br.readLine()) != null) {

					z[i] = new complex(Float.parseFloat(line), 0);

					i++;

				}// end while

				if (f.equals("pulseData.txt")) {

					System.out.println("Working");
					for (int j = i; j < z.length; j++) {

						z[j] = new complex(0, 0);

					}// end j for

				}// end if

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return z;

	}

	/**
	 * Prints the two d compelx array
	 * 
	 * @param a
	 *            the two d complex array
	 */
	public void printtwod(complex[][] a) {

		complex[][] z = new complex[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				z[i][j] = a[i][j];

			}// end j for

		}// end i for

		for (int i = 0; i < z.length; i++) {

			for (int j = 0; j < z[0].length; j++) {

				System.out.print("[" + z[i][j].x + " :  " + z[i][j].y + "]");
				// System.out.println("[" + b[i][j].x + " :  " + b[i][j].y +
				// "]");

			}// end j for

			System.out.println();

		}// end i for

	}// end printtwod
}// end class
