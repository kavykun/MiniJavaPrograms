import java.awt.geom.Point2D;

public class Points2DLabel {

	char label;
	Point2D point;
	int index;
	String neighbor;

	public Points2DLabel(Point2D p, char l, int i) {

		point = p;
		label = l;
		index = i;

	}// end constructor Points

	public char getLabel() {
		return label;
	}

	public void setLabel(char label) {
		this.label = label;
	}

	public Point2D getPoint() {
		return point;
	}

	public void setPoint(Point2D point) {
		this.point = point;
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	public String getNeighbor() {
		return neighbor;
	}

	public void setNeighbor(String neighbor) {
		this.neighbor = neighbor;
	}

}// end Points
