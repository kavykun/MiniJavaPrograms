public class complex {

	double x, y;
	ErrorHandler error;
	complex z, u;

	public complex(float d, float b) {

		this.x = d;
		this.y = b;

	}
	public complex cnum(float x, float y) {

		u = new complex(0, 0);
		u.x = x;
		u.y = y;
		return u;
	}

	public complex cadd(complex u, complex v) {

		z = new complex(0, 0);
		z.x = u.x + v.x;
		z.y = u.y + v.y;
		return z;
	}

	public complex csub(complex u, complex v) {

		z = new complex(0, 0);
		z.x = u.x - v.x;
		z.y = u.y - v.y;
		return z;
	}

	public complex cmul(complex u, complex v) {

		z = new complex(0, 0);
		z.x = u.x * v.x - u.y * v.y;
		z.y = u.x * v.y + u.y * v.x;
		return z;
	}

	public complex cdiv(complex u, complex v) {

		z = new complex(0, 0);
		float d = (float) (v.x * v.x + v.y * v.y), eps = epsilon();
		if (d < eps * eps)
			error.nerror(104, "cdiv", 0, 0, "");
		z.x = (u.x * v.x + u.y * v.y) / d;
		z.y = (u.y * v.x - u.x * v.y) / d;
		return z;
	}


	public complex csqrt(complex u) {

		z = new complex(0, 0);
		float r = (float) Math.sqrt(cmag(u)), phi = cphi(u) / 2;
		double x, y;

		x = r * Math.cos(phi);
		z.x = (float) x;
		y = r * Math.sin(phi);
		z.y = (float) y;
		return z;
	}

	public void negateY() {

		this.y *= -1;

	}

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

}// end class
