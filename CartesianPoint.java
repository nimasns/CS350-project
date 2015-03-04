
public class CartesianPoint implements Comparable<CartesianPoint>{
	int x = 0;
	int y = 0;
	
	CartesianPoint(int x_coord, int y_coord) {
		this.x = x_coord;
		this.y = y_coord;
	}

	@Override
	public int compareTo(CartesianPoint pt) {
		if (this.x < pt.x) return -1;
		else if (this.x == pt.x && this.y < pt.y) return -1;
		else if (this.x == pt.x && this.y == pt.y) return 0;
		else if (this.x == pt.x && this.y > pt.y) return 1; 
		else return 1;
	}
}
