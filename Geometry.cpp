#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

class Figure;
class Point;
class Vector;
class Line;
class Segment;
class Ray;
class Polygon;



class Figure {
public:
	virtual void shift(const Vector&) = 0;
	virtual bool includes(const Point&) const = 0;
	virtual bool cross(const Segment&) const = 0;
};



class Point: Figure {
public:
	double x_;
	double y_;

	double get_x() const { return x_; }
	double get_y() const { return y_; }

	void set_x(double val) { x_ = val; }		 
	void set_y(double val) { y_ = val; }		 

	Point(double x = 0.0, double y = 0.0) : x_(x), y_(y) {}
	void shift(const Vector& vec);
	double dist(const Point& value) const;
	bool includes(const Point&) const;
	bool cross(const Segment&) const;
};



class Vector {
public:
	double x_;
	double y_;

	Vector(double x = 0.0, double y = 0.0) : x_(x), y_(y) {}
	Vector(const Point& value) : x_(value.x_), y_(value.y_) {}
	Vector(const Vector& value) : x_(value.x_) , y_(value.y_) {}
	Vector(const Point& a, const Point& b);

	double get_x() const { return x_; }
	double get_y() const { return y_; }
	double modul() const { return sqrt(y_ * y_ + x_ * x_); }
	void shift(const Vector& vec) {
		this->x_ += vec.x_;
		this->y_ += vec.y_;
	}
	double operator* (const Vector& vec) const{
		return x_ * vec.x_ + y_ * vec.y_;
	}
	double operator^ (const Vector& vec) const{
		return x_ * vec.y_ - y_ * vec.x_;
	}
	double triangle(const Vector& vec) {
		return (*this ^ vec) / 2;
	}
	Vector operator+ (const Vector& vec) const {
		return Vector(x_ + vec.x_, y_ + vec.y_);
	}
	Vector normal() {
		return Vector(-y_, x_);
	}
	double angle(const Vector&) const;
};		



class Line : public Figure {	
public:
	double a_;
	double b_;
	double c_;

	Line(double a, double b, double c) : a_(a), b_(b), c_(c) {}
	Line(const Point& a, const Point& b);
	Line(const Point& a, const Vector& vec);


	double get_a() const { return a_; }
	double get_b() const { return b_; }
	double get_c() const { return c_; }

	Vector get_vec() const;
	bool bool_cross(const Line&) const;
	double dist_line(const Line&) const;
	Point point_cross(const Line&) const;
	Point get_point() const;
	Point projection(const Point&) const;
	double dist_point(const Point&) const;

	void shift(const Vector&);
	bool includes(const Point&) const;
	bool cross(const Segment&) const;
};


class Segment : public Figure {
public:
	Point a_;
	Point b_;

	Segment(const Point& a, const Point& b) : a_(a), b_(b) {}
	void shift(const Vector&);
	bool includes(const Point&) const;
	bool cross(const Segment&) const;

	Point get_a() const { return a_;}
	Point get_b() const { return b_;}
	Line get_line() const;
	double distance_to_point(const Point&) const;
	double distance_to_segment(const Segment&) const;
};



class Ray : public Figure {
public:
	Point a_;
	Vector vec_;

	Ray(const Point& a, const Vector& vec) : a_(a), vec_(vec) {}
	Ray(const Point&, const Point&);

	Point get_a() const { return a_;}
	Vector get_vec() const { return vec_;}
	void shift(const Vector&);
	bool includes(const Point&) const;
	bool cross(const Segment&) const;
	Line get_line() const { return Line(a_, vec_);}
	double distance_to_point(const Point&) const;
};


class Polygon : public Figure {
public:
	int n_;
	Point* values_;

	Polygon(int n, Point *mas) {
		n_ = n;
		values_ = mas;
	}
	
	void shift(const Vector&);
	bool includes(const Point&) const;
	bool cross(const Segment&) const;
	bool prominent() const;
	double square() const;
	~Polygon() { delete [] values_; }
};



int main() {
}



double Point::dist(const Point& value) const {
	Vector tmp(*this, value);
	return tmp.modul();
}

void Point::shift(const Vector& vec) {
	x_ += vec.x_;
	y_ += vec.y_;
}

bool Point::includes(const Point& value) const {
	if ((x_ == value.x_) && (y_ == value.y_))
		return true;
	else
		return false;
}

bool Point::cross(const Segment& seg) const {
	return seg.includes(*this);
}

double Vector::angle(const Vector& vec) const {
	return acos(*this * vec / (this->modul() * vec.modul()));
}

Vector::Vector(const Point& a, const Point& b) : x_(b.x_ - a.x_), y_(b.y_ - a.y_) {}



Line::Line(const Point& a, const Point& b) {
	if (a.x_ == b.x_) {
		a_ = 1; b_ = 0; c_ = -b.x_;
	}
	else {
		if (a.y_ == b.y_) {
			a_ = 0; b_ = 1; c_ = -b.y_;
		}
		else {
			a_ = 1 / (a.x_ - b.x_);
			b_ = 1 / (b.y_ - a.y_);
			c_ = a.y_ / (a.y_ - b.y_) - a.x_ / (a.x_ - b.x_);
		}
	}
}

Line::Line(const Point& a, const Vector& vec) {
		Point b = a;
		b.set_x(b.x_ + vec.x_);
		b.set_y(b.y_ + vec.y_);
		
		if (a.x_ == b.x_) {
			a_ = 1; b_ = 0; c_ = -b.x_;
		}
		else {
			if (a.y_ == b.y_) {
				a_ = 0; b_ = 1; c_ = -b.y_;
			}
			else {
				a_ = 1 / (a.x_ - b.x_);
				b_ = 1 / (b.y_ - a.y_);
				c_ = a.y_ / (a.y_ - b.y_) - a.x_ / (a.x_ - b.x_);
			}
		}
	}

Point Line::get_point() const {
	if (a_ == 0)
		return Point(1, -c_ / b_);
	if (b_ == 0)
		return Point(-c_ / a_, 1);
	return Point(1, (-c_ - a_) / b_);
}

void Line::shift(const Vector& vec) {
	c_ = c_ - vec.x_ * a_ - vec.y_ * b_;
	return;
}

bool Line::includes(const Point& value) const {
	double tmp = a_ * value.x_ + b_ * value.y_;
	return ((tmp <= -c_ + 0.0000000000001) && (tmp >= -c_ - 0.0000000000001));
}

bool Line::cross(const Segment & value) const {
	Vector vec = this->get_vec();
	Point pnt = this->get_point();
	Vector a(pnt, value.a_);
	Vector b(pnt, value.b_);
	return (this->includes(value.a_) || this->includes(value.b_) || (vec ^ a) * (vec ^ b) <= 0);
}

bool Line::bool_cross(const Line& that) const {
	if (this->get_vec() ^ that.get_vec()) 
		return true;
	else return 0;
}

Point Line::point_cross(const Line& that) const {
	double x, y;	
	if (b_ == 0) {
		x = -c_ / a_;
		y = (-that.c_ - that.a_ * x) / that.b_;
		Point ans(x, y);
		return ans;
	}
	if (that.b_ == 0) {
		x = -that.c_ / that.a_;
		y = (-c_ - a_ * x) / b_;
		Point ans(x, y);
		return ans;
	}
	x = -(c_ * that.b_ - that.c_ * b_) / (that.b_ * a_ - b_ * that.a_);
	y = (-c_ - a_ * x) / b_;
	Point ans(x, y);
	return ans;
}

double Line::dist_line(const Line& that) const {
	if ( ((a_ * that.b_) == (b_ * that.a_)) && ((c_ * that.b_) == (b_ * that.c_)) && ((a_ * that.c_) == (c_ * that.a_)) )
		return 0;
	if (a_ == 0)
		return fabs(c_ / b_ - that.c_ / that.b_);
	if (b_ == 0)
		return fabs(c_ / a_ - that.c_ / that.a_);
	return fabs(c_ - that.c_ * a_ / that.a_) / sqrt(a_ * a_ + b_ * b_);
}

Vector Line::get_vec() const {
	return Vector(-b_, a_);
}

Point Line::projection(const Point& value) const {
	return this->point_cross(Line(value, get_vec().normal()));
}

double Line::dist_point(const Point& value) const {
	return value.dist(projection(value));
}



Line Segment::get_line() const {
	Line l(this->a_, this->b_);
	return l;
}

void Segment::shift(const Vector& vec) {
	a_.shift(vec);
	b_.shift(vec);
}

bool Segment::includes(const Point& value) const {
	Ray r1(a_, b_);
	Ray r2(b_, a_);
	return ((r1.includes(value)) && (r2.includes(value)));
}

bool Segment::cross(const Segment& seg) const {
	Line l1(a_, b_), l2(seg.a_, seg.b_);
	double x1 = min(a_.x_, b_.x_);
	double x2 = max(a_.x_, b_.x_);
	double x3 = min(seg.a_.x_, seg.b_.x_);
	double x4 = max(seg.a_.x_, seg.b_.x_);
	double y1 = min(a_.y_, b_.y_);
	double y2 = max(a_.y_, b_.y_);
	double y3 = min(seg.a_.y_, seg.b_.y_);
	double y4 = max(seg.a_.y_, seg.b_.y_);
	bool tmp = ((x1 <= x4) && (x2 >= x3) && (y1 <= y4) && (y2 >= y3));
	return (l1.cross(seg) && l2.cross(*this) && tmp);
}

double Segment::distance_to_point(const Point& value) const {
	if (this->includes(this->get_line().projection(value)))
		return get_line().dist_point(value);
	else
		return min(value.dist(a_), value.dist(b_));
}

double Segment::distance_to_segment(const Segment& seg) const {
	return min(min(distance_to_point(seg.a_), distance_to_point(seg.b_)),
		min(seg.distance_to_point(this->a_), seg.distance_to_point(this->b_)));
}



Ray::Ray(const Point& a, const Point& b) {
	a_ = a;
	vec_ = Vector(a, b);
}

void Ray::shift(const Vector& vec) {
	a_.shift(vec);
	return;
}

bool Ray::includes(const Point& value) const {
	return (this->get_line().includes(value) && ((value.x_ - this->a_.x_) * this->vec_.x_ >= 0) && ((value.y_ - this->a_.y_) * this->vec_.y_ >= 0)) ;
}

bool Ray::cross(const Segment& seg) const {
	return ((this->get_line().cross(seg)) && (this->includes(this->get_line().bool_cross(seg.get_line()))));
}

double Ray::distance_to_point(const Point& value) const {
	if (this->includes(this->get_line().projection(value)))
		return get_line().dist_point(value);
	else
		return a_.dist(value);
}



void Polygon::shift(const Vector& vec) {
	for (int i = 0; i < n_; ++i) {
		values_[i].shift(vec);
	}
	return;
}

bool Polygon::includes(const Point& value) const{
	if (Segment(values_[n_ - 1], values_[0]).includes(value)) {
		return true;
	}
	for (int i = 0; i < n_ - 1; i++) {
		if (Segment(values_[i], values_[i + 1]).includes(value)) {
			return true;
		}
	}	
	double angles = 0;
	Vector one(value, values_[n_ - 1]);
	Vector two(value, values_[0]);
	int k = (one ^ two) < 0 ? 1 : -1;
	angles += (one.angle(two)  * k);
	for (int i = 0; i < n_ - 1; i++) {
		one = two;
		two = Vector(values_[i + 1].x_ - value.x_, values_[i + 1].y_ - value.y_);
		k = (one ^ two) < 0 ? 1 : -1;
		angles += (one.angle(two) * k);
	}
	return !(fabs(angles) < 1.0000000001);
}

bool Polygon::cross(const Segment& seg) const {
	for (int i = 0; i < n_ - 1; i++) {
		if (Segment(values_[i], values_[i + 1]).cross(seg))
			return true;
	}
	if (Segment(values_[0], values_[n_ - 1]).cross(seg))
		return true;
	if (this->includes(seg.a_) || this->includes(seg.b_) )
		return true;
	return false;
}

bool Polygon::prominent() const{
	Vector one(values_[n_ - 1], values_[0]);
	Vector two(values_[0], values_[1]);

	double mult = one ^ two;
	for (int i = 0; i < n_ - 1; i++) {
		if (i == n_ - 2) {
			one = Vector(values_[n_ - 2], values_[n_ - 1]);
			two = Vector(values_[n_ - 1], values_[0]);
		}
		else {
			one = Vector(values_[i], values_[i + 1]);
			two = Vector(values_[i + 1], values_[i + 2]);
		}
		if (mult == 0) {
			mult = one ^ two;
		}
		if ((one ^ two) * mult < 0) {
			return false;
		}
	}
	return true;
}

double Polygon::square() const {
	double square = 0;
	square += Vector(values_[n_ - 1]).triangle(Vector(values_[n_ - 1], values_[0]));
	for (int i = 0; i < n_ - 1; ++i) {
		square += Vector(values_[i]).triangle(Vector(values_[i], values_[i + 1]));
	}
	return fabs(square);
}



double max(double a, double b) {
	return a > b ? a : b;
}

double min(double a, double b) {
	return a < b ? a : b;
}