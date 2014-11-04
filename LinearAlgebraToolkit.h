#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

#ifndef LINEARALGEBRATOOLKIT_H
#define LINEARALGEBRATOOLKIT_H

//the concept of a vector is going to be handled fairly loosely here in that 
//I'm allowing them to switch freely between a row and a column vector. You 
//can treat this as essentially automatically transposing a vector to make it 
//work with whatever operation it needs to be used in at the time. This means
//that I'm allowing vector*matrix and matrix*vector with the same vector, regardless
//of the fact that this doesn't make sense from a purely mathematical standpoint. 
//Just pretend I'm doing a transpose operation and live happily in your mathematical bliss.
class Vector {

	private:
		int size;
		double* vals;

	public:
		Vector() : size(0), vals(NULL) {}
		Vector(int s) : size(s) { 
				allocateMemory();
				fillVector(0);
		}
		//2d case
		Vector(double a, double b) : size(2) {
			allocateMemory();
			vals[0] = a; vals[1] = b;
		}
		//3d case
		Vector(double a, double b, double c) : size(3) {
			allocateMemory()
			vals[0] = a; vals[1] = b; vals[2] = c;
		}
		//I assume you're giving the correct size of the 
		//array along with the array. It will break if you don't. 
		Vector(int s, double* arr) : size(s) {
				allocateMemory();
				for(int i = 0; i < s; i++) 
					vals[i] = arr[i];
			}
		}
		
		//copy constructor
		Vector(const Vector& v) {
			if(size != v.size) {
				destroy();
				size = v.size;
				allocateMemory();
			}
			for(int i = 0; i < size; i++) 
				vals[i] = v.vals[i];
		}

		Vector(const Vector* v) {
			if(size != v->size) {
				destroy();
				size = v->size;
				allocateMemory();
			}
			for(int i = 0; i < size; i++) 
				vals[i] = v->vals[i];
		}
		
		//destructor
		~Vector() {
			destroy();
		}
		
        //callable destructor
        void destroy() {
			if(vals != NULL) {
				delete [] vals;
			}
			vals = NULL;
			size = 0;
		}

        void allocateMemory() {
            if(size <= 0) {
                size = 0;
                vals = NULL;
            } else {
                vals = new double[size];
            }
        }

		int getSize() {
			return size;
		}

		//sets the value at pos
		bool setVal(int pos, double val) {
			if(pos < 0 || pos >= size) 
				return false;
			vals[pos] = val;
			return true;
		}

        //overwrites the values in the vector with the values from the supplied vector
		void setVals(Vector v) {
			if(size != v.size) {
				destroy();
				size = v.size;
			}
			allocateMemory();
			for(int i = 0; i < v.size; i++) {
				vals[i] = v[i];
			}
		}

		//overwrites the values in the vector using the values from the supplied array
		void setVals(int size, double* arr) {
			setVals(Vector(size, arr));
		}
		
		bool is2d() {
			if(size == 2) 
				return true;
			return false;
		}
		
		bool is3d() {
			if(size == 3) 
				return true;
			return false;
		}

		double getVal(int index) {
			if(index < 0 || index >= size) {
				return -12345;
			}
			return vals[index];
		}
		
        double operator[](int index) {
			return getVal(index);
		}

		//scalar multiplication
		//Note: I included functionality for both Vector * scalar
		//		   and scalar * Vector. they may or may not make
		//		   sense mathematically, but intuitively, it makes 
		//		   sense to me, so....I'm leaving it.
		friend Vector operator*(Vector v, double scale) {
			Vector result = new Vector(v.size, v.vals);
			for(int i = 0; i < v.size; i++) {
				result.vals[i] *= scale;
			}
			return result;
		}
		
		friend Vector operator*(double scale, Vector v) {
			return v*scale;
		}

		Vector operator*=(double scale) {
			for(int i = 0; i < size; i++) {
				vals[i] *= scale;
			}
			return *this;
		}
		
		//scalar division
		//Note: I only included Vector/scalar because it doesn't
		//		  make much sense to divide a number by a vector
		friend Vector operator/(Vector v, double scale) {
			return v*(1.0/scale);
		}
		
		Vector operator/=(double scale) {
			for(int i = 0; i < size; i++) {
				vals[i] /= scale;
			}
			return *this;
		}

		//dot product
		friend double operator*(Vector u, Vector v) {
			if(u.size != v.size) 
				return -12345;
			double sum = 0;
			for(int i = 0; i < u.size; i++) {
				sum += u.vals[i] * v.vals[i];
			}
			return sum;
		}

		friend double operator*(Vector u, Vector* v) {
			return u*(*v);
		}
		
		friend double operator*(Vector* u, Vector v) {
			return (*u)*v;
		}
		
		//more dot product in case people don't know that * is the dot product
		static double dot(Vector u, Vector v) {
			return u*v;
		}
		
		static double dot(Vector* u, Vector v) {
			return u*v;
		}
		
		static double dot(Vector u, Vector* v) {
			return u*v;
		}
		
		static double dot(Vector* u, Vector* v) {
			return (*u)*(*v);
		}
		
        //so much dot product...
		double dot(Vector v) {
			return (*this)*v;
		}
		
		double dot(Vector* v) {
			return (*this)*(*v);
		}

		//cross product. Vectors must be 3 dimensional to get
		//a meaningful result. 
		static Vector cross(Vector u, Vector v) {
			if(!(u.is3d() && v.is3d())) 
				return Vector();
			Vector result(3);
			//formula for cross product of two vectors;
			result.vals[0] = u[1]*v[2] - u[2]*v[1];
			result.vals[1] = u[2]*v[0] - u[0]*v[2];
			result.vals[2] = u[0]*v[1] - u[1]*v[0];
			return result;
		}
		
		static Vector cross(Vector* u, Vector v) {
			return cross(*u, v);
		}
		
		static Vector cross(Vector* u, Vector* v) {
			return cross(*u, *v);
		}

		static Vector cross(Vector u, Vector* v) {
			return cross(u, *v);
		}
		
		Vector cross(Vector v) {
			return cross(*this, v);
		}

		Vector cross(Vector* v) {
			return cross(*this, *v);
		}

		//vector addition
		friend Vector operator+(Vector u, Vector v) {
			if(u.size != v.size) 
				return Vector();
			Vector result(u.size);
			for(int i = 0; i < u.size; i++) {
				result.vals[i] = u.vals[i] + v.vals[i];
			}
			return result;
		}

		friend Vector operator+(Vector u, Vector* v) {
			return u+(*v);
		}

		friend Vector operator+(Vector* u, Vector v) {
			return (*u)+v;
		}

        Vector operator+=(Vector v) {
            *this = *this + v;
            return *this;
        }

		//vector subtraction. You'll note this is just adding
		//the negative version of the vector on the right. 
		friend Vector operator-(Vector u, Vector v) {
			return (u + (v*-1));
		}
		
		friend Vector operator-(Vector u, Vector* v) {
			return (u + ((*v)*-1));
		}
	
		friend Vector operator-(Vector* u, Vector v) {
			return ((*u) + (v*-1));
		}

        Vector operator-=(Vector v) {
            *this = *this - v;
            return *this;
        }
	
		//equals operator
		Vector& operator=(Vector v) {
			if(this == &v) {
				return *this;
			}
            destroy();
			size = v.size;
			allocateMemory();
			for(int i = 0; i < size; i++) { 
				vals[i] = v[i];
			}
			return *this;
		}

		Vector& operator=(Vector* v) {
			if(this == v) {
				return *this;
			}
            destroy();
			size = v->size;
            allocateMemory();
			for(int i = 0; i < size; i++) { 
				vals[i] = v->vals[i];
			}
			return *this;
		}

        //comparators. 
		bool operator==(Vector v) {
			if(v.size != size) {
				return false;
			}
			for(int i = 0; i < size; i++) {
				if(v[i] != vals[i]) {
					return false;
				}
			}
			return true;
		}
		
		bool operator!=(Vector v) {
			return !(*this == v);
		}
		
        //this is where a less than and greater than comparator would go,
        //but I don't think that makes a lot of sense for vectors. you could
        //be comparing magnitude, position, angle, or any number of things,
        //so if you want to compare those, you'll have to do it yourself

		friend ostream& operator<<(ostream& out, Vector v) {
			if(v.size > 0) {
				out << "[ " << v[0];
				for(int i = 1; i < v.size; i++) {
					out << ", " << v[i];
				}
				out << " ]";
			} else {
				out << "[ ]";
			}
			return out;
		}
		
        //expected format: [ 1.0, 2, 4.3 ]
        //also, this overwrites the vector with the vector being entered
		friend istream& operator>>(istream& in, Vector& v) {
			v.destroy();
			char garbage;
			double val;
			in >> garbage;
			while(garbage != ']') {
				in >> val >> garbage;
				v.pushBack(val);
			}
			return in;
		}

		//returns the length of the vector
		double length() {
			double sum = 0;
			for(int i = 0; i < size; i++) {
				sum += (vals[i] * vals[i]);
			}
			return sqrt(sum);
		}
		
		static double length(Vector v) {
			return v.length();
		}
		
		//same thing as length...I've heard it called magnitude before though too.
		double magnitude() {
			return length();
		}
		
		static double magnitude(Vector v) {
			return v.length();
		}
		
		//normalize the vector this is called on
		void normalize() {
			*this = (*this)/length();
		}
		
		//returns a normalized version of the vector this was called on
		Vector getNormalizedVector() {
			Vector result(size, vals);
			return result.normalize();
		}
		
		//I've heard it called the unit vector more than normalized vector, so this is for that...
		void toUnitVector() {
           normalize();
        }

        Vector getUnitVector() {
			return getNormalizedVector();
		}
		
        //returns the largest number in the vector
		double max() {
			if(size == 0) 
				return 0;
			double max = vals[0];
			for(int i = 1; i < size; i++) {
				if(vals[i] > max)
					max = vals[i];
			}
			return max;
		}
		
        //returns the smallest number in the vector
		double min() {
			if(size == 0) 
				return 0;
			double min = vals[0];
			for(int i = 1; i < size; i++) {
				if(vals[i] < min)
					min = vals[i];
			}
			return min;
		}
		
		bool isUnit() {
			if(length() == 1) {
				return true;
			}
			return false;
		}

		void clear() {
			destroy();
		}

        //changes the size of the vector to the newSize, keeping as many of the values as 
        //possible in the case of a smaller size, and adding zeroes in the case of a 
        //larger size
		void changeSize(int newSize) {
			if(newSize <= 0) {
                destroy();
                return;
            }
            double* temp = vals;
			vals = new double[newSize];
			if(newSize < size) {
				for(int i = 0; i < newSize; i++) {
					vals[i] = temp[i];
				}
			} else {
				for(int i = 0; i < size; i++) {
					vals[i] = temp[i];
				}
				for(int i = size; i < newSize; i++) {
					vals[i] = 0;
				}
			}
			size = newSize;
			delete [] temp;
		}

		void pushBack(double newVal) {
			changeSize(size+1);
			setVal(size-1, newVal);
		}

		static Vector zeroVector(int size) {
			return Vector(size);
		}

		static void fillVector(Vector v, double fillNumber) {
			v.fillVector(fillNumber);
		}
		
        //fills the vector with the given fill number
		void fillVector(double fillNumber) {
			for(int i = 0; i < size; i++) {
				vals[i] = fillNumber;
			}
		}

        //returns true if there are any values that are not zero in the vector
		bool isNonZero() {
			for(int i = 0; i < size; i++) {
				if(vals[i] != 0) {
					return true;
				}
			}
			return false;
		}

        //returns true if all values in the vector are zero, or if the size is zero
		bool isZero() {
			return !isNonZero();
		}

        //gives the index of the first value that isn't zero, or in the case 
		//of all values being zero, returns the size. mainly added this to 
        //make reducing matrices to row-echelon form easier
		int getFirstNonZeroElementIndex() {
			for(int i = 0; i < size; i++) {
				if(vals[i] != 0) {
					return i;
				}
			}
			return size;
		}

        //returns an array with the same values as the vector
		double* toArray() {
			double* result = new double[size];
			for(int i = 0; i < size; i++) {
				result[i] = vals[i];
			}
			return result;
		}

        //this is meant for applying a rotation matrix to a vector,
        //so I'm assuming it's a vector 3 when this is called, and 
        //that the last value should be a 1. it will work for any
        //size vector, but it won't be particularly useful.
        Vector toVector4() {
            changeSize(4);
            setVal(3, 1.0);
            return *this;
        }

        Vector toVector3() {
            changeSize(3);
            return *this;
        }

        bool isEmpty() {
            if(size == 0) {
                return true;
            }
            return false;
        }
		
		friend Vector abs(Vector v) {
			Vector result;
			for(int i = 0; i < v.size; i++)  {
				result.pushBack(abs(v[i]));
			}
			return result;
		}
};

//pretty much all the comments made for the vector stuff hold for the point stuff
class Point {

	private: 
		int size;
		double* vals;
	
	public: 
	
		Point() : size(0), vals(NULL) {}
		Point(int s) : size(s) { allocateMemory(); zeroPoint(); }
		Point(double a, double b) : size(2) { allocateMemory(); vals[0] = a; vals[1] =b; }
		Point(double a, double b,  double c) : size(3) {
			allocateMemory();
			vals[0] = a;
			vals[1] = b;
			vals[2] = c;
		}
		Point(int s, double* arr) : size(s) {
			allocateMemory();
			for(int i = 0; i < size; i++) {
				vals[i] = arr[i];
			}
		}
		
		Point(const Point& p) {
			size = p.size;
			vals = new double[size];
			for(int i = 0; i < size; i++) {
				vals[i] = p.vals[i];
			}
		}

		Point(const Point* p) {
			size = p->size;
			vals = new double[size];
			for(int i = 0; i < size; i++) {
				vals[i] = p->vals[i];
			}
		}
		
		~Point() {
			deallocateMemory();
		}
		
        void destroy() {
			deallocateMemory();
		}

        void allocateMemory() {
			if(size <= 0) {
                size = 0; 
                vals = NULL;
            } else {
                vals = new double[size];
            }
        }
		
		void deallocateMemory() {
			if(vals != NULL) {
				delete [] vals;
			}
			vals = NULL;
			size = 0;
		}

		int getSize() {
			return size;
		}
		
		double operator[](int index) {
			return getVal(index);
		}
		
		double getVal(int index) {
			if(index < 0 || index > size) {
				return -12345;
			}
			return vals[index];
		}
		
		void setVal(int index, double val) {
			if(index >= 0 && index < size) {
				vals[index] = val;
			}
		}
		
        void setVals(Point p) {
            if(size != p.size) {
                destroy();
				size = p.size;
				allocateMemory();
			}
			for(int i = 0; i < size; i++) {
				vals[i] = p[i];
			}
        }

		void setVals(int size, double* arr) {
			setVals(Point(size, arr));
		}
		
		bool is2d() {
			if(size == 2) 
				return true;
			return false;
		}
		
		bool is3d() {
			if(size == 3) 
				return true;
			return false;
		}
		
		Point& operator=(Point p) {
			if(this == &p) {
				return *this;
			}
			destroy();
			size = p.size;
            allocateMemory();
            for(int i = 0; i < size; i++) { 
				vals[i] = p.vals[i];
			}
			return this;
		}

		Point& operator=(Point* p) {
			if(this == p) {
				return *this;
			}
			destroy();
			size = p->size;
            allocateMemory();
            for(int i = 0; i < size; i++) { 
				vals[i] = p->vals[i];
			}
			return *this;
		}

		bool operator==(Point p) {
			if(p.size != size) {
				return false;
			}
			for(int i = 0; i < size; i++) {
				if(p[i] != vals[i]) {
					return false;
				}
			}
			return true;
		}
		
		bool operator!=(Point p) {
			return !(*this == p);
		}
		
		friend Point operator*(Point p, double s) {
			Point result(p.size, p.vals);
			for(int i = 0; i < p.size; i++) {
				result.vals[i] *= s;
			}
			return result;
		}
		
		friend Point operator*(double s, Point p) {
			return p*s;
		}

		Point operator*=(double s) {
			for(int i = 0; i < size; i++) {
				vals[i] *= s;
			}
			return *this;
		}
		
		friend Point operator/(double s, Point p) {
			return p*(1.0/s);
		}
		
		friend Point operator/(Point p, double s) {
			return p*(1.0/s);
		}
		
		Point operator/=(double s) {
			for(int i = 0; i < size; i++) {
				vals[i] /= s;
			}
			return *this;
		}
		
		friend Point operator+(Point p, Vector v) {
			if(p.size != v.getSize()) {
				return Point();
			}
			Point result(p.size);
			for(int i = 0; i < p.size; i++) {
				result.vals[i] = p.vals[i] + v.getVal(i);
			}
			return result;
		}
		
		friend Point operator+(Point p, Vector* v) {
			return p+(*v);
		}

        Point operator+=(Vector v) {
            *this = *this + v;
            return *this;
        }
		
		friend Point operator-(Point p, Vector v) {
			return p + (v*-1);
		}
		
		friend Point operator-(Point p, Vector* v) {
			return p + ((*v)*-1);
		}

        Point operator-=(Vector v) {
            *this = *this - v;
            return *this;
        }
		
		friend Vector operator-(Point p, Point q) {
			if(p.size != q.size) {
				return Vector();
			}
			Vector result(p.size);
			for(int i = 0; i < p.size; i++) {
				result.setVal(i, p[i] - q[i]);
			}
			return result;
		}
		
		friend Vector operator-(Point* p, Point q) {
			return (*p)-q;
		}
		
		friend Vector operator-(Point p, Point* q) {
			return p-(*q);
		}
		
		friend ostream& operator<<(ostream& out, Point p) {
			if(p.size > 0) {
				out << "[ " << p[0];
				for(int i = 1; i < p.size; i++) {
					out << ", " << p[i];
				}
				out << " ]";
			} else {
				out << " [ ] ";
			}
			return out;
		}
		
		friend istream& operator>>(istream& in, Point& p) {
			p.destroy();
			char garbage;
			double val;
			in >> garbage;
			while(garbage != ']') {
				in >> val >> garbage;
				p.pushBack(val);
			}
			return in;
		}
		
		void zeroPoint() {
			for(int i = 0; i < size; i++) {
				vals[i] = 0;
			}
		}
		
		void fillPoint(double fillNumber) {
			for(int i = 0; i < size; i++) {
				vals[i] = fillNumber;
			}
		}
		
		void clear() {
			zeroPoint();
		}
		
		double max() {
			if(size == 0) 
				return 0;
			double max = vals[0];
			for(int i = 1; i < size; i++) {
				if(vals[i] > max)
					max = vals[i];
			}
			return max;
		}
		
		double min() {
			if(size == 0) 
				return 0;
			double min = vals[0];
			for(int i = 1; i < size; i++) {
				if(vals[i] < min)
					min = vals[i];
			}
			return min;
		}
		
		void changeSize(int newSize) {
			double* temp = vals;
			vals = new double[newSize];
			if(newSize < size) {
				for(int i = 0; i < newSize; i++) {
					vals[i] = temp[i];
				}
			} else {
				for(int i = 0; i < size; i++) {
					vals[i] = temp[i];
				}
				for(int i = size; i < newSize; i++) {
					vals[i] = 0;
				}
			}
			size = newSize;
			delete [] temp;
		}
	
		void pushBack(double newVal) {
			changeSize(size+1);
			setVal(size-1, newVal);
		}
		
		double* toArray() {
			double* result = new double[size];
			for(int i = 0; i < size; i++) {
				result[i] = vals[i];
			}
			return result;
		}

        Vector toVector() {
            return Vector(size, vals);
        }

        bool isEmpty() {
            if(size == 0) {
                return true;
            }
            return false;
        }
		
		friend Point abs(Point p) {
			Point result;	
			for(int i = 0; i < p.size; i++) {
				result.pushBack(abs(v[i]));
			}
			return result;
		}
};
/********************MATRIX************************/
class Matrix {

	private:
		int rows, cols;
		double** mat;
	
	public:
	
		Matrix() : rows(0), cols(0), mat(NULL) {}
		Matrix(int r, int c) : rows(r), cols(c) {
			allocateMemory();
			zeroMatrix();
		}

		Matrix(int r, int c, double init) : rows(r), cols(c) {
			allocateMemory();
			fillMatrix(init);
		}

		Matrix(const Matrix& m) {
			rows = m.rows;
			cols = m.cols;
			allocateMemory();
			for(int i = 0; i < rows; i++) {
				for(int j = 0; j < cols; j++) {
					mat[i][j] = m.mat[i][j];
				}
			}
		}

		Matrix(const Matrix* m) {
			rows = m->rows;
			cols = m->cols;
			allocateMemory();
			for(int i = 0; i < rows; i++) {
				for(int j = 0; j < cols; j++) {
					mat[i][j] = m->mat[i][j];
				}
			}
		}

		~Matrix() { deallocateMemory(); }

		Matrix& operator=(Matrix m) {
			if(this == &m) {
				return *this;
			}
			deallocateMemory();
			rows = m.rows;
			cols = m.cols;
			allocateMemory();
			for(int i = 0; i < rows; i++) {
				for(int j = 0; j < cols; j++) {
					mat[i][j] = m.mat[i][j];
				}
			}
		}
		
		Matrix& operator=(Matrix* m) {
			if(this == m) {
				return *this;
			}
			deallocateMemory();
			rows = m->rows;
			cols = m->cols;
			allocateMemory();
			for(int i = 0; i < rows; i++) {
				for(int j = 0; j < cols; j++) {
					mat[i][j] = m->mat[i][j];
				}
			}
		}

		void allocateMemory() {
            if(rows == 0 || cols == 0) {
                mat = NULL;
                rows = 0;
                cols = 0;
                return;
            }
			mat = new double*[rows];
			for(int i = 0; i < rows; i++) {
				mat[i] = new double[cols];
			}
		}

		void deallocateMemory() {
			if(mat != NULL) {
				for(int i = 0; i < rows; i++) {
					delete [] mat[i];
				}
				delete [] mat;
			}
		}
		
		void clear() {
			deallocateMemory();
		}

		double operator()(int r, int c) {
			if(!isInBounds(r,c)) { 
				return -12345;
			} else {
				return mat[r][c];
			}
		}

		//if I return a vector, it lets you access items in a m[i][j] fashion...loopholes!
		Vector operator[](int r) {
			return Vector(cols, mat[r]);
		}

		int numRows() { return rows; }
		int numcols() { return cols; }
		int rowLength() { return cols; }
		int colLength() { return rows; }

		void zeroMatrix() {
			fillMatrix(0);
		}

		void fillMatrix(double fill) {
			for(int i = 0; i < rows; i++) {
				for(int j = 0; j < cols; j++) {
					mat[i][j] = fill;
				}
			}
		}

		bool isInBounds(int r, int c) {
			if(r >= rows || c >= cols || r < 0 || c < 0) {
				return false;
			} else {
				return true;
			}
		}
		
		void setVal(int r, int c, double val) {
			if(isInBounds(r,c)) {
				mat[r][c] = val;
			}
		}

		void setRow(int rowNum, double val) {
			if(isInBounds(rowNum, 0)) {
				for(int i = 0; i < cols; i++) {
					mat[rowNum][i] = val;
				}
			}
		}

		void setRow(int rowNum, Vector v) {
			if(isInBounds(rowNum, 0)) {
				for(int i = 0; i < cols; i++) {
					mat[rowNum][i] = v[i];
				}
			}
		}

		void setCol(int colNum, double val) {
			if(isInBounds(0, colNum)) {
				for(int i = 0; i < rows; i++) {
					mat[i][colNum] = val;
				}
			}
		}
		
		void setCol(int colNum, Vector v) {
			if(isInBounds(0, colNum)) {
				for(int i = 0; i < rows; i++) {
					mat[i][colNum] = v[i];
				}
			}
		}

		bool isSquare() {
			if(rows == cols) 
				return true;
			return false;
		}

        //this adds a row to the matrix. it assumes that if the matrix has 
        //no values in it, that the matrix was meant to be empty, and thus makes
        //a 1xN matrix with the given vector as the only row. if you'd like to keep
        //all the rows of zeroes, see addRowKeepingZeroes()
		bool addRow(Vector v) {
			if(!isZeroMatrix()) {
				if(v.getSize() != cols) {
					return false;
				}
				Matrix temp = *this;
				deallocateMemory();
				rows++;
				allocateMemory();
				for(int i = 0; i < rows-1; i++) {
					for(int j = 0; j < cols; j++) {
						mat[i][j] = temp.mat[i][j];
					}
				}
			} else {
				deallocateMemory();
				cols = v.getSize();
				rows = 1;
				allocateMemory();
			}
			for(int i = 0; i < cols; i++) {
				mat[rows-1][i] = v[i];
			}
			return true;
		}

		bool addRow(int size, double* arr) {
			return addRow(Vector(size, arr));
		}

        bool addRowKeepingZeroes(Vector v) {
            if(!isZeroMatrix()) {
                return addRow(v);
            } else {
                if(v.getSize() != rowLength()) {
                    return false;
                }
                int tempRows = rows;
                int tempCols = cols;
                deallocateMemory();
                rows = tempRows+1;
                cols = tempCols;
                allocateMemory();
                zeroMatrix();
                setRow(rows-1, v);
                return true;
            }
        }

        bool addRowKeepingZeroes(int s, double* arr) {
            return addRowKeepingZeroes(Vector(s, arr));
        }
		
		bool addCol(Vector v) {
			if(!isZeroMatrix()) {
				if(v.getSize() != rows) {
					return false;
				}
				Matrix temp = this;
				deallocateMemory();
				cols++;
				allocateMemory();
				for(int i = 0; i < rows; i++) {
					for(int j = 0; j < cols-1; j++) {
						mat[i][j] = temp.mat[i][j];
					}
					mat[i][cols-1] = v[i];
				}
			} else {
				deallocateMemory();
				rows = v.getSize();
				cols = 1;
				allocateMemory();
				for(int i = 0; i < rows; i++) {
					mat[i][0] = v[i];
				}
			}
			return true;
		}

		bool addCol(int size, double* arr) {
			return addCol(Vector(size, arr));
		}

        bool addColKeepingZeroes(Vector v) {
            if(!isZeroMatrix()) {
                return addCol(v);
            } else {
                if(v.getSize() != colLength()) {
                    return false;
                }
                int tempRows = rows;
                int tempCols = cols;
                deallocateMemory();
                rows = tempRows;
                cols = tempCols+1;
                allocateMemory();
                zeroMatrix();
                setCol(cols-1, v);
                return true;
            }
        }

        bool addColKeepingZeroes(int s, double* arr) {
            return addColKeepingZeroes(Vector(s, arr));
        }

        //returns true if the dimensions of the input matrix are the same as the
        //dimensions of the calling object
		bool compareDimensions(Matrix m) {
			return (m.rows == rows && m.cols == cols) 
		}

		bool compareDimensions(Matrix* m) {
			return (m->rows == rows && m->cols == cols) 
		}

		static bool compareDimensions(Matrix m, Matrix l) {
			return (m.rows == l.rows && m.cols == l.cols)
		}

		bool deleteRow(int rowNum) {
			if(!isInBounds(rowNum, 0)) {
				return false;
			}
			Matrix temp = this;
			deallocateMemory();
			rows--;
			allocateMemory();
			bool skippedRow = false;
			for(int i = 0; i < temp.rows; i++) {
				if(i == rowNum) {
					skippedRow = true;
					continue;
				}
				for(int j = 0; j < cols; j++) {
					if(skippedRow)
						mat[i-1][j] = temp.mat[i][j];
					else 
						mat[i][j] = temp.mat[i][j];
				}
			}
			return true;
		}

		bool deleteCol(int colNum) {
			if(!isInBounds(0, colNum)) {
				return false;
			}
			Matrix temp = this;
			deallocateMemory();
			cols--;
			allocateMemory();
			for(int i = 0; i < rows; i++) {
				bool skippedCol = false;
				for(int j = 0; j < temp.cols; j++) {
					if(j == colNum) {
						skippedCol = true;
						continue;
					}
					if(skippedCol)
						mat[i][j-1] = temp.mat[i][j];
					else
						mat[i][j] = temp.mat[i][j];
				}
			}
			return true;
		}

		friend ostream& operator<<(ostream& out, Matrix m) {
			out.setf(ios::fixed);
			out.precision(3);
			out << "[ ";
			for(int i = 0; i < m.colLength(); i++) {
				if(i != 0)
					out << "  ";
				out << Vector(m.rowLength(), m.mat[i]) << endl;
			}
			out << "  ]" << endl;
			return out;
		}
		
		friend istream& operator>>(istream& in, Matrix& m) {
			m.clear();
			char garbage;
			Vector v;
			in >> garbage;
			while(garbage != ']') {
				in >> v >> garbage;
				m.addRow(v);
			}
			return in;
		}

		friend Matrix operator+(Matrix m, Matrix l) {
			if(!m.compareDimensions(l))
				return Matrix();
			Matrix result(m.rows, m.cols);
			for(int i = 0; i < result.rows; i++) {
				for(int j = 0; j < result.cols; j++) {
					result.mat[i][j] = m.mat[i][j] + l.mat[i][j];
				}
			}
			return result;
		}

		friend Matrix operator+(Matrix* m, Matrix l) {
			return (*m)+l;
		}
		 
		friend Matrix operator+(Matrix m, Matrix* l) {
			return m+(*l);
		}

		friend Matrix operator-(Matrix m, Matrix l) {
			if(!m.compareDimensions(l))
				return Matrix();
			Matrix result(m.rows, m.cols);
			for(int i = 0; i < result.rows; i++) {
				for(int j = 0; j < result.cols; j++) {
					result.mat[i][j] = m.mat[i][j] - l.mat[i][j];
				}
			}
			return result;
		}

		friend Matrix operator-(Matrix* m, Matrix l) {
			return (*m)-l;
		}
		 
		friend Matrix operator-(Matrix m, Matrix* l) {
			return m-(*l);
		}

		friend Matrix operator*(Matrix m, Matrix l) {
			if(m.rowLength() != l.colLength())
				return Matrix();
			Matrix result(m.colLength(), l.rowLength(), 0);
			for(int i = 0; i < result.rows; i++) {
				for(int j = 0; j < result.cols; j++) {
					for(int k = 0; k < m.cols; k++) {
						result.mat[i][j] += m.mat[i][k]*l.mat[k][j];
					}
				}
			}
			return result;
		}

		friend Matrix operator*(Matrix* m, Matrix l) {
			return (*m)*l;
		}

		friend Matrix operator*(Matrix m, Matrix* l) {
			return m*(*l);
		}

        //for a row vector
        friend Vector operator*(Vector v, Matrix m) {
            if(v.getSize() != m.colLength()) {
                return Vector();
            }
            Vector result(m.rowLength());
            for(int i = 0; i < m.rowLength(); i++) {
                double sum = 0.0;
                for(int j = 0; j < m.colLength(); j++) {
                    sum += v[j]*m[j][i];
                }
                result.setVal(i, sum);
            }
            return result;
        }

        friend Vector operator*(Vector* v, Matrix m) {
            return (*v)*m;
        }


        //for a column vector
		friend Vector operator*(Matrix m, Vector v) {
			if(v.getSize() != m.rowLength()) {
                return Vector();
            }
            Vector result(m.colLength());
            for(int i = 0; i < m.colLength(); i++) {
                double sum = 0.0;
                for(int j = 0; j < m.rowLength(); j++) {
                    sum += v[j]*m[i][j];
                }
                result.setVal(i, sum);
            }
            return result;
		}

		friend Vector operator*(Matrix m, Vector* v) {
			return m*(*v);
		}

		friend Point operator*(Point p, Matrix m) {
			if(p.getSize() != m.colLength()) {
				return Point();
			}
			Point result(m.rowLength());
			for(int i = 0; i < m.rowLength(); i++) {
				double sum = 0.0;
				for(int j = 0; j < m.colLength(); j++) {
					sum += p[j]*m[j][i];
				}
				result.setVal(i, sum);
			}
			return result;
		}

		friend Point operator*(Point* p, Matrix m) {
			return (*p)*m;
		}
		
		friend Point operator*(Matrix m, Point p) {
			if(p.getSize() != m.rowLength()) {
                return Point();
            }
            Point result(m.colLength());
            for(int i = 0; i < m.colLength(); i++) {
                double sum = 0.0;
                for(int j = 0; j < m.rowLength(); j++) {
                    sum += p[j]*m[i][j];
                }
                result.setVal(i, sum);
            }
            return result;
		}

		friend Point operator*(Matrix m, Point* p) {
			return m*(*p);
		}
		
		friend Matrix operator*(Matrix m, double s) {
			Matrix result(m.rows, m.cols);
			for(int i = 0; i < result.rows; i++) {
				for(int j = 0; j < result.cols; j++) {
					result.mat[i][j] = m.mat[i][j]*s;
				}
			}
			return result;
		}

		friend Matrix operator/(Matrix m, double s) {
			return m*(1.0/s);
		}

		bool operator==(Matrix m) {
			for(int i = 0; i < rows; i++) {
				for(int j = 0; j < cols; j++) {
                    if(mat[i][j] != m[i][j]) {
                        return false;
					}
				}
			}
			return true;
		}

		bool operator!=(Matrix m) {
			return !(*this == m);
		}

		Matrix getRowEchelonForm() {
			Matrix result;
			result = *this;
			if(result.colLength() <= 1) {
				return result;
			}
			result.orderRows();

			for(int i = 0; i < result.colLength(); i++) {
				Vector currRow = result.getRowAsVector(i);
				int counter = currRow.getFirstNonZeroElementIndex();                
				if(counter >= result.rowLength()) {
					continue;
				}
				for(int j = i+1; j < result.colLength(); j++) {
					Vector temp = result.getRowAsVector(j);
					double numer = temp.getVal(counter);
					if(numer == 0) {
						continue;
					} else {
						double denom = currRow[counter];
						double mult = numer/denom;
						temp = temp - currRow*mult;
						result.setRow(j, temp);
					}
				}
				result.orderRows();
			}

			return result;
		}

		Matrix getReducedRowEchelonForm() {
			Matrix result = getRowEchelonForm();

			for(int i = 0; i < result.colLength(); i++) {
				Vector temp = result.getRowAsVector(i);
				int ind = temp.getFirstNonZeroElementIndex();
				temp /= temp[ind];
				result.setRow(i, temp);
			}

			for(int i = 1; i < result.colLength(); i++) {
				Vector currRow = result.getRowAsVector(i);
				int ind = currRow.getFirstNonZeroElementIndex();
				for(int j = i-1; j >= 0; j--) {
					Vector temp = result.getRowAsVector(j);
					temp = temp - currRow*temp[ind];
					result.setRow(j, temp);
				}
			}
			return result;
		}

		//this function orders the rows: 
		//                              all zero rows go to the bottom
		//                              rows with leading zeroes go below rows without
		void orderRows() {
			if(colLength() <= 1 || isOrdered()) {
				return;
			}
			//get all the zero rows to the bottom
			int numZeroRows = 0;
			for(int i = 0; i < colLength() - numZeroRows; i++) {
				Vector temp = getRowAsVector(i);
				while(temp.isZero()) {
					swapRows(i, colLength() - 1 - numZeroRows);
					numZeroRows++;
					temp = getRowAsVector(i);
				}
			}
			//order the remaining rows...this is terribly inefficient
			//should probably figure out a good way to do this instead of just throwing it
			//in a while !ordered loop...yay bubble sort. 
			while(!isOrdered()) {
				for(int i = 0; i < colLength() - numZeroRows; i++) {
					Vector temp = getRowAsVector(i);
					int ind = temp.getFirstNonZeroElementIndex();
					for(int j = i+1; j < colLength() - numZeroRows; j++) {
						Vector temp2 = getRowAsVector(j);
						int ind2 = temp2.getFirstNonZeroElementIndex();
						if(ind2 < ind) {
							swapRows(i, j);
							break;
						}
					}
				}
			}

		}

		Matrix reduceToRowEchelonForm() {
			*this = getRowEchelonForm();
			return *this;
		}

		Matrix reduceToReducedRowEchelonForm() {
			*this = getReducedRowEchelonForm();
			return *this;
		}

		bool isInRowEchelonForm() {
			if(colLength() <= 1) {
				return false;
			}
			//checking all non-zero rows are above any zero rows
			int firstZeroRow = -1;
			for(int i = 0; i < colLength(); i++) {
				Vector temp = getRowAsVector(i);
				if(temp.isZero()) {
					firstZeroRow = i;
					break;
                }
			}
			int lastNonZeroRow = -1;
			for(int i = 0; i < colLength(); i++) {
				Vector temp = getRowAsVector(i);
				if(temp.isNonZero()) {
					lastNonZeroRow = i;
				}
			}
			if(firstZeroRow != -1 && firstZeroRow < lastNonZeroRow) {
				return false;
			}
			//checking all leading coefficients are to the right of the leading coefficient of the row above it
			int prevCol = -1;
			int currCol = 0;
			for(int i = 0; i < colLength(); i++) {
				Vector temp = getRowAsVector(i);
				if(temp.isNonZero()) {
					for(int j = 0; j < temp.getSize(); j++) {
						if(temp[j] != 0) {
							currCol = j;
							if(currCol <= prevCol) {
								return false;
							}
							prevCol = currCol;
							break;
						}
					}
                }
			}
			return true;
		}

		bool isInReducedRowEchelonForm() {
			if(!isInRowEchelonForm()) {
				return false;
			}
			//checking if the leading coefficient is one and that the leading coefficient is the only nonzero entry in the column
			for(int i = 0; i < colLength(); i++) {
				Vector tempRow = getRowAsVector(i);
				if(tempRow.isZero()) {
					continue;
				}
				int leadingCoefficient = -1; 
				for(int j = 0; j < tempRow.getSize(); j++) {
					if(tempRow[j] != 0) {
						if(tempRow[j] != 1) {
							return false;
						}
						leadingCoefficient = j;
						break;
					}
				}
				Vector tempCol = getColAsVector(leadingCoefficient);
				for(int j = 0; j < tempCol.getSize(); j++) {
					if(tempCol.getVal(j) != 0 && j != i) {
						return false;
					}
				}
			}
			return true;
		}

		bool isOrdered() {
			for(int i = 0; i < colLength(); i++) {
				Vector temp = getRowAsVector(i);
				int ind = temp.getFirstNonZeroElementIndex();
				for(int j = i+1; j < colLength(); j++) {
					Vector temp2 = getRowAsVector(j);
					int ind2 = temp2.getFirstNonZeroElementIndex();
					if(ind2 < ind) {
						return false;
					}
				}
			}
			return true;
		}

		Matrix inverse() {
			if(!isSquare() || determinant() == 0) {
				return Matrix();
			}
			Matrix result(rows, cols);
			result = *this;
            //add identity matrix to the right side of the matrix
			for(int i = 0; i < rowLength(); i++) {
				Vector v;
				for(int j = 0; j < colLength(); j++) {
				   if(i == j) {
					   v.pushBack(1);
				   } else {
					   v.pushBack(0);
				   }
				}
				result.addCol(v);
			}

            //the reduced row echelon form of the matrix if the determinant is not zero
            //and the matrix is square will be the identity matrix, leaving the inverse on
            //the right side of the new matrix we made
			result.reduceToReducedRowEchelonForm();

            //get rid of the identity matrix on the left side
			for(int i = rowLength() - 1; i >= 0; i--) {
				result.deleteCol(i);
			}

			return result;
		}

		Matrix inversify() {
			*this = inverse();
			return *this;
		}

		//run time for this is n!/2...only use if you have to and store the value when you can as to not re-run it.
		//TODO: make this not shitty run time.
		double determinant() {
			if(!isSquare()) {
				return -12345;
			}
			if(rows == 2 && cols == 2) {
				return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
			}
			double sum = 0;
			double sign = 1;
			for(int i = 0; i < rowLength(); i++) {
				Matrix temp = createSubMatrix((*this), i);
				sum += sign*mat[0][i]*temp.determinant();
				sign*=-1;
			}
			return sum;
		}
		
		//helper function for determinant calculations
		private Matrix createSubMatrix(Matrix mat, int pos) {
			Matrix temp(mat.colLength()-1, mat.rowLength()-1);
			bool skipped = false;
			for(int j = 0; j < mat.rowLength(); j++) {
				Vector v;
				if(j != pos) {
					for(int k = 1; k < mat.colLength(); k++) {
						v.pushBack(mat[k][j]);
					}
				} else {
					skipped = true;
				}
				if(v.getSize() != 0) {
					if(skipped) {
						temp.setCol(j-1, v);
					} else {
						temp.setCol(j, v);
					}
				}
			}
			return temp;
		}

		Vector getRowAsVector(int rowNum) {
			if(!isInBounds(rowNum, 0)) {
				return Vector();
			}
			Vector result;
			for(int i = 0; i < rowLength(); i++) {
				result.pushBack(mat[rowNum][i]);
			}
			return result;
		}

		double* getRowAsArray(int rowNum) {
			return getRowAsVector(rowNum).toArray();
		}

		Vector getColAsVector(int colNum) {
			if(!isInBounds(0, colNum)) {
				return Vector();
			}
			Vector result;
			for(int i = 0; i < colLength(); i++) {
				result.pushBack(mat[i][colNum]);
			}
			return result;
		}

		double* getColAsArray(int colNum) {
			return getColAsVector(colNum).toArray();
		}

		bool swapRows(int row1, int row2) {
			if(!isInBounds(row1, 0) || !isInBounds(row2, 0)) {
				return false;
			}
			
			Vector temp1 = getRowAsVector(row1);
			Vector temp2 = getRowAsVector(row2);

			setRow(row1, temp2);
			setRow(row2, temp1);

			return true;
		}

		bool swapCols(int col1, int col2) {
			if(!isInBounds(0, col1) || !isInBounds(0, col2)) {
				return false;
			}

			Vector temp1 = getColAsVector(col1);
			Vector temp2 = getColAsVector(col2);

			setCol(col1, temp2);
			setCol(col2, temp1);

			return true;
		}

		//get the matrix starting at startRow, and ending at endCol. both parameters are inclusive
		Matrix getPartialMatrixByCols(int startCol, int endCol) {
			if(!isInBounds(0, startCol) || !isInBounds(0, endCol)) {
				return Matrix();
			}
			Matrix result;

			for(int i = startCol; i <= endCol; i++) {
				Vector temp = getColAsVector(i);
				result.addCol(temp);
			}

			return result;
		}

		Matrix getPartialMatrixByRows(int startRow, int endRow) {
			if(!isInBounds(startRow, 0) || !isInBounds(endRow, 0)) {
				return Matrix();
			}
			Matrix result;

			for(int i = startRow; i <= endRow; i++) {
				Vector temp = getRowAsVector(i);
				result.addRow(temp);
			}

			return result;
		}

		bool deleteRows(int startRow, int endRow) {
			if(!isInBounds(startRow, 0) || !isInBounds(endRow, 0)) {
				return false;
			}
			//could optimize this to copy the rows left over if we're 
			//deleting more than half the matrix. 
			for(int i = startRow; i <= endRow; i++) {
				deleteRow(i);
			}
		}

		bool deleteCols(int startCol, int endCol) {
			if(!isInBounds(0, startCol) || !isInBounds(0, endCol)) {
				return false;
			}
			//could optimize this to copy the columns left over if we're 
			//deleting more than half the matrix. 
			for(int i = startCol; i <= endCol; i++) {
				deleteCol(i);
			}
		}

		bool isZeroMatrix() {
			for(int i = 0; i < rows; i++) {
				for(int j = 0; j < cols; j++) {
					if(mat[i][j] != 0) {
						return false;
					}
				}
			}
			return true;
		}

		//solves the Ax=b equation, assuming the calling object is A, and b is provided.
		//returns x;
		Vector solveAxb(Vector b) {
			if(!isSquare() || determinant() == 0) {
				return Vector();
			}
			Vector x(rowLength());
			x = inverse() * b;

			return x;
		}

		Matrix getTranspose() {
			Matrix result(cols, rows);
			for(int i = 0; i < rows; i++) {
				for(int j = 0; j < cols; j++) {
					result.mat[j][i] = mat[i][j];
				}
			}
			return result;
		}

		Matrix transpose() {
			*this = getTranspose();
			return *this;
		}

		static Matrix getIdentityMatrix(int size) {
			Matrix result(size, size);
			for(int i = 0; i < size; i++) {
				result.mat[i][i] = 1;
			}
			return result;
		}

		Vector toVector() {
			if(rows != 1 || cols != 1) {
				return Vector();
			}
			if(rows == 1) {
				return getColAsVector(0);    
			} else {
				return getRowAsVector(0);
			}
		}

		bool isSingular() {
			if(determinant() == 0) {
				return true; 
			} else {
				return false;
			}
		}

        bool applyRotationInRadians(double rads, string axis) {
            if(rows == 2 && cols == 2) {
                Matrix rot;
                rot.addRow(Vector(cos(rads), -1*sin(rads)));
                rot.addRow(Vector(sin(rads), cos(rads)));
                *this = rot*(*this);
                return true;
            } else if(rows == 3 && cols == 3) {
                Matrix rot;
                if(axis == "x") {
                    rot.addRow(Vector(1,0,0));
                    rot.addRow(Vector(0,cos(rads), -1*sin(rads)));
                    rot.addRow(Vector(0,sin(rads), cos(rads)));
                } else if (axis == "y") {
                    rot.addRow(Vector(cos(rads), 0, sin(rads)));
                    rot.addRow(Vector(0, 1, 0));
                    rot.addRow(Vector(-1*sin(rads), 0, cos(rads)));
                } else if (axis == "z") {
                    rot.addRow(Vector(cos(rads), -1*sin(rads), 0));
                    rot.addRow(Vector(sin(rads), cos(rads), 0));
                    rot.addRow(Vector(0,0,1));
                }
                *this = rot*(*this);
                return true;
            }
            return false;
        }

        bool applyRotationInDegrees(double degrees, string axis) {
            double radians = degrees*3.14159/180.0;
            return applyRotationInRadians(radians, axis);
        }

        //trans should be in the form of [ x, y ] for translating a 2x2
        //and in the form of [ x, y, z ] for translating a 3x3
        bool applyTranslation(Vector trans) {
            if(rows == 2 && cols == 2) {
                trans.pushBack(1);
                Matrix ident = getIdentityMatrix(3);
                ident.setCol(2, trans);
                changeSize(rows+1, cols+1);
                setRow(rows, Vector(1,1,1));
                *this = ident*(*this);
                deleteRow(rows);
                deleteCol(cols);
                return true;
            } else if (rows == 3 && cols == 3) {
                trans.pushBack(1);
                Matrix ident = getIdentityMatrix(4);
                ident.setCol(3, trans);
                changeSize(rows+1, cols+1);
                Vector v(4);
                v.fillVector(1);
                setRow(rows, v);
                *this = ident*(*this);
                deleteRow(rows-1);
                deleteCol(cols-1);
                return true;
            }
            return false;
        }

        void changeSize(int numRows, int numCols) {
            while(rows > numRows) {
                deleteRow(rows-1);
            }
            while(cols > numCols) {
                deleteCol(cols-1);
            }
            while(rows < numRows) {
                addRow(Vector(cols));
            }
            while(cols < numCols) {
                addCol(Vector(rows));
            }
        }

        Matrix get2dRotationMatrix(double rads) {
            Matrix rot;
            rot.addRow(Vector(cos(rads), -1*sin(rads)));
            rot.addRow(Vector(sin(rads), cos(rads)));
            return rot;
        }
        
        //returns a matrix that will rotate another matrix 'rads' redians about the 'axis' axis.
        //this is designed purely for rotation, however, and does not provide the capability of
        //adding a translation component without some work on your part (read: adding another row 
        //and column to the result of this function call)
        Matrix get3dRotationMatrix(double rads, string axis) {
            Matrix rot;
            if(axis == "x") {
                rot.addRow(Vector(1,0,0));
                rot.addRow(Vector(0,cos(rads), -1*sin(rads)));
                rot.addRow(Vector(0,sin(rads), cos(rads)));
            } else if (axis == "y") {
                rot.addRow(Vector(cos(rads), 0, sin(rads)));
                rot.addRow(Vector(0, 1, 0));
                rot.addRow(Vector(-1*sin(rads), 0, cos(rads)));
            } else if (axis == "z") {
                rot.addRow(Vector(cos(rads), -1*sin(rads), 0));
                rot.addRow(Vector(sin(rads), cos(rads), 0));
                rot.addRow(Vector(0,0,1));
            }
            return rot;
        }

        bool isEmpty() {
            if(rows == 0 && cols == 0) {
                return true;
            }
            return false;
        }

        friend Matrix abs(Matrix m) {
            Matrix result = m;
            for(int i = 0; i < result.rows; i++) {
                for(int j = 0; j < result.cols; j++) {
                    result.mat[i][j] = abs(result.mat[i][j]);
                }
            }
            return result;
        }
};

#endif //LINEARALGEBRATOOLKIT_H