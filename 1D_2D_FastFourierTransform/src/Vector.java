/**
 * Class to keep each vector
 * 
 * @author Kavy
 * 
 */
public class Vector {

	double x, y;
	public double X = 0;
	public double Y = 0;
	public char label;

	public Vector(double a, double b) {

		this.x = a;
		this.y = b;

	}// end constructor vector

	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}

	public Vector subtract(Vector a, Vector b) {

		Vector result;
		double F = (a.getX() - b.getX());
		double L = (a.getY() - b.getY());
		result = new Vector(F, L);

		return result;

	}

	/**
	 * A to string method for a vector.
	 * 
	 * @return
	 */
	public String tostring() {
		String ex = Double.toString(X);
		String why = Double.toString(Y);
		String ans = (label + " " + "[" + x + "] [" + y + "]");
		return ans;

	}

}// end class vector
