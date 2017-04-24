#include<iostream>
#include<cstring>
#include<cstdlib>
#include<stdexcept>
#include <cmath>

using namespace std;


//Rational


int nod(int a, int b){
    return b ? nod(b, a % b) : a;
}

class RationalDivisionByZero : std::logic_error {
public:
    RationalDivisionByZero();
};

RationalDivisionByZero::RationalDivisionByZero() : logic_error("error") {}


class Rational {
private:
    int p;
    int q;

    void reduce() {
        int a = nod(this->p, this->q);
        this->p /= a;
        this->q /= a;
        if (this->q < 0){
            this->p *= -1;
            this->q *= -1;
        }
    }

public:
    Rational(): p(0), q(1) {}
    Rational(int p, int q): p(p), q(q) {
        this->reduce();
    }
    Rational(int a): p(a), q(1) {}
    Rational(const Rational& other): p(other.p), q(other.q) {}

    int getNumerator() const {
        return this->p;
    }
    int getDenominator() const {
        return this->q;
    }

	Rational operator+() const{
        return Rational(this->p, this->q);
    }
    Rational operator+(const Rational& rat) const;
	Rational& operator+=(const Rational& rat){
		*this = *this + rat;
        return *this;
    }

    Rational operator-() const{
        return Rational(-this->p, this->q);
    }
    Rational operator-(const Rational& rat) const {
		return *this + (-rat);
    }
    Rational operator-(const int a) const {
		return *this - Rational(a);    
    }
	Rational& operator-=(const Rational& rat){
		*this = *this - rat;
        return *this;
    }

    Rational operator*(const Rational& rat) const;
    Rational operator*(const int a) const {
		return *this * Rational(a);  
    }
	Rational& operator*=(const Rational& rat){
        *this = *this * rat;
        return *this;
    }

    Rational operator/(const Rational& rat) const;
    Rational operator/(const int a) const{
		return *this / Rational(a);
    }    
    Rational& operator/=(const Rational& rat){
        *this = *this / rat;
        return *this;
    }

    bool operator==(const Rational& rat) const{
        return (this->p == rat.p && this->q == rat.q);
    }
    bool operator!=(const Rational& rat) const{
        return !(*this == rat);
    }
    bool operator>(const Rational& rat) const{
        return (*this - rat).p > 0;
    }
    bool operator<(const Rational& rat) const{
        return (*this - rat).p < 0;
    }
    bool operator>=(const Rational& rat) const{
        return (*this - rat).p >= 0;
    }
    bool operator<=(const Rational& rat) const{
        return (*this - rat).p <= 0;
    }

    Rational& operator++(){
        this->p += this->q;
        return *this;
    }
    Rational operator++(int){
        Rational old(*this);
        ++(*this);
        return old;
    }

    Rational& operator--(){
        this->p -= this->q;
        return *this;
    }
    Rational operator--(int not_used){
        Rational old(*this);
        --(*this);
        return old;
    }
   
    friend Rational operator+(const int, const Rational&);
    friend Rational operator-(const int, const Rational&);
    friend Rational operator*(const int, const Rational&);
    friend Rational operator/(const int, const Rational&);
    friend bool operator>(const int, const Rational&);
    friend bool operator<(const int, const Rational&);
    friend bool operator>=(const int, const Rational&);
    friend bool operator<=(const int, const Rational&);
    friend bool operator==(const int, const Rational&);
    friend bool operator!=(const int, const Rational&);
    friend istream& operator>>(istream&, Rational&);
	friend ostream& operator<<(ostream&, const Rational&);
};



Rational Rational::operator/(const Rational& rat) const {
    if (rat.p == 0)
        throw RationalDivisionByZero();
    else {
        int nod1 = nod(this->p, rat.p), nod2 = nod(this->q, rat.q);
        Rational res(this->p / nod1 * rat.q / nod2, this->q / nod2 * rat.p / nod1);
        return res;
    }
}

Rational Rational::operator*(const Rational& rat) const {
	int nod1 = nod(this->p, rat.q), nod2 = nod(this->q, rat.p);
	Rational res(this->p / nod1 * rat.p / nod2, this->q / nod2 * rat.q / nod1);
    return res;
}

Rational Rational::operator+(const Rational& rat) const {
        int nod1 = nod(this->q, rat.q);
        int denom = this->q / nod1 * rat.q;
        int nom = rat.q / nod1 * this->p + this->q / nod1 * rat.p;
        Rational res(nom, denom);
        res.reduce();
        return res;
    }

Rational operator+(const int a, const Rational& rat) {
    return Rational(a) + rat;
}

Rational operator-(const int a, const Rational& rat) {
    return Rational(a) - rat;
}

Rational operator*(const int a, const Rational& rat) {
    return Rational(a) * rat;
}

Rational operator/(const int a, const Rational& rat) {
    return Rational(a) / rat;
}

bool operator>(const int a, const Rational& rat){
    return rat < a;
}

bool operator<(const int a, const Rational& rat){
    return rat > a;
}

bool operator>=(const int a, const Rational& rat){
    return rat <= a;
}

bool operator<=(const int a, const Rational& rat){
    return rat >= a;
}

bool operator==(const int a, const Rational& rat){
    return rat == a;
}

bool operator!=(const int a, const Rational& rat){
    return rat != a;
}

istream& operator>>(std::istream& is, Rational& rat) {
	char str[50];
	is >> str;
    char* split = strtok(str, "/");
    rat.p = atoi(split);
    split = strtok(NULL, "/");
	rat.q = split ? atoi(split) : 1;
    rat.reduce();
    return is;
}

ostream& operator<<(ostream &os, const Rational& rat) {
    os << rat.p;
    if (rat.q != 1)
        os << "/" << rat.q;
    return os;
}


// Matrix


class MatrixWrongSizeError : logic_error {
public:
    MatrixWrongSizeError();
};

MatrixWrongSizeError::MatrixWrongSizeError() : logic_error("Error") {}

class MatrixIndexError : logic_error {
public:
    MatrixIndexError();
};

MatrixIndexError::MatrixIndexError() : logic_error("Error") {}

class MatrixIsDegenerateError : logic_error {
public:
    MatrixIsDegenerateError();
};

MatrixIsDegenerateError::MatrixIsDegenerateError() : logic_error("Error") {}



template <typename T>
class Matrix {
private:
    int row_count;
    int col_count;
    T **values;
protected:
    Matrix();
    void assign(const Matrix& that);
    void set_mem(const int col_count, const int raw_count);
    void delete_mem();
    void swap_rows(const int, const int);
    void row_add(const int, const int, const T);
public:
    Matrix(const Matrix<T>& );
    Matrix(int, int);
    ~Matrix();
    T get(int, int) const;
    void set(int, int, T);   
    int get_rows_numb() const;
    int get_columns_numb()  const;
    Matrix operator-(const Matrix& ) const;
    Matrix operator+(const Matrix& ) const;
    Matrix operator*(const Matrix& ) const;
    Matrix& operator=(const Matrix& );
    Matrix& operator+=(const Matrix& );
    Matrix& operator-=(const Matrix& );
    Matrix& operator*=(const Matrix& );
    Matrix& operator*=(const T k);
    T operator()(int, int) const;
    Matrix operator*(const T k) const;
    Matrix operator/(const T k) const;
    Matrix get_trans () const;
    Matrix& transpose ();
    
	template <typename C>
    friend istream& operator>>(istream& in, const Matrix<C>& a);
    template <typename C>
    friend ostream& operator<<(ostream& res, const Matrix<C>& a);
    template <typename C>
    friend Matrix<C> operator*(const C k, const Matrix<C>& a);
};


template <typename T>
void Matrix<T>::delete_mem() {
    for (int i = 0; i < this->row_count; i++)
        delete []this->values[i];
    delete []this->values;
}
template <typename T>
void Matrix<T>::set_mem(const int row_count, const int col_count) {
    this->row_count = row_count;
    this->col_count = col_count;
    this->values = new T *[this->row_count];
    for (int i = 0; i < this->row_count; i++)
    {
        this->values[i] = new T[this->col_count];
    }
}
template <typename T>
void Matrix<T>::assign(const Matrix<T>& that){
    set_mem(that.row_count, that.col_count);
    for (int i = 0; i < this->row_count; i++){
        for (int j = 0; j < this->col_count; j++){
            this->values[i][j] = that.values[i][j];
        }
    }
}
template <typename T>
Matrix<T>::Matrix(){}
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& that){
    assign(that);
}
template <typename T>
Matrix<T>::Matrix(int row_count, int col_count){
    set_mem(row_count, col_count);
    for (int i = 0; i < row_count; i++){
        for (int j = 0; j < col_count; j++){
            this->values[i][j] = 0;
        }
    }
}
template <typename T>
Matrix<T>::~Matrix() {
    delete_mem();
}
template <typename T>
void Matrix<T>::set(int i, int j, T val){
    this->values[i][j] = val;
}
template <typename T>
T Matrix<T>::get(int i, int j) const{
    return this->values[i][j];
}
template <typename T>
T Matrix<T>::operator()(int row, int column) const{
    if (row >= this->row_count || column >= this->col_count) {
        throw MatrixIndexError();
    }
    return this->values[row][column];
}
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& that) const{
    if (this->row_count != that.row_count || this->col_count != that.col_count) {
        throw MatrixWrongSizeError();
    }
    Matrix res(this->row_count, this->col_count);
    for (int i = 0; i < this->row_count; i++){
        for (int j = 0; j < this->col_count; j++){
            res.values[i][j] = this->values[i][j] + that.values[i][j];
        }
    }
    return res;
}
template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& that) const {
    if (this->row_count != that.row_count || this->col_count != that.col_count) {
        throw MatrixWrongSizeError();
    }
    return *this + that * (-1);
}
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& that) const {
    if (this->col_count != that.row_count) {
        throw MatrixWrongSizeError();
    }
    Matrix res(this->row_count, that.col_count);
    for (int i = 0; i < res.row_count; i++) {
        for (int j = 0; j < res.col_count; j++) {
            res.values[i][j] = 0;
            for (int q = 0; q < this->col_count; q++) {
                res.values[i][j] += this->values[i][q] * that.values[q][j];
            }
        }
    }
    return res;
}
template <typename T>
Matrix<T> Matrix<T>::operator*(const T k) const {
    Matrix res(this->row_count, this->col_count);
    for (int i = 0; i < this->row_count; i++) {
        for (int j = 0; j < this->col_count; j++) {
            res.values[i][j] = this->values[i][j] * k;
        }
    }
    return res;
}
template <typename T>
Matrix<T> Matrix<T>::operator/(const T k) const {
    Matrix res(this->row_count, this->col_count);
    for (int i = 0; i < this->row_count; i++) {
        for (int j = 0; j < this->col_count; j++) {
            res.values[i][j] = this->values[i][j] / k;
        }
    }
    return res;
}
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& that) {
    if (this != &that) {
        delete_mem();
        assign(that);
    }
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix& that) {
    *this = *this + that;
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix& that) {
    *this = *this - that;
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix& that) {
    *this = *this * that;
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const T k) {
    *this = *this * k;
    return *this;
}
template <typename T>
Matrix<T> Matrix<T>::get_trans () const {
    Matrix res(this->col_count, this->row_count);
    for (int i = 0; i < this->row_count; i++) {
        for (int j = 0; j < this->col_count; j++) {
            res.values[j][i] = this->values[i][j];
        }
    }
    return res;
}
template <typename T>
Matrix<T>& Matrix<T>::transpose() {
    *this = get_trans();
    return *this;
}

template <typename T>
int Matrix<T>::get_rows_numb() const{
    return this->row_count;
}
template <typename T>
int Matrix<T>::get_columns_numb()  const {
    return  this->col_count;
}

template <typename T>
void Matrix<T>::swap_rows(const int first, const int second) {
    for (int i = 0; i < this->col_count; i++) {
        T temp = this->values[first][i];
        this->values[first][i] = this->values[second][i];
        this->values[second][i] = -temp;
    }
}

template <typename T>
void Matrix<T>::row_add(const int first_row, const int second_row, const T k) {
    for (int i = 0; i < this->col_count; i++) {
        this->set(first_row, i, this->values[first_row][i] + this->values[second_row][i] * k);
    }
}
template <typename C>
istream& operator>>(istream &in, const Matrix<C>& a){
    for (int i = 0; i < a.row_count; i++){
        for (int j = 0; j < a.col_count; j++)
            in >> a.values[i][j];
    }
    return in;
}
template <typename C>
ostream& operator<<(ostream& res, const Matrix<C>& a){
    for (int i = 0; i < a.row_count; i++){
        for (int j = 0; j < a.col_count; j++)
            res << a.values[i][j] << " ";
        res << endl;
    }
    return res;
}
template <typename C>
Matrix<C> operator*(const C k, const Matrix<C>& a) {
    return a * k;
}

//Square matrix

template <typename T>
class SquareMatrix : public Matrix<T> {
private:
    void assign(const Matrix<T>& that);
    SquareMatrix get_minor(const int i, const int j) const;
public:
    SquareMatrix (int size);
    SquareMatrix ();
    SquareMatrix (const Matrix<T>& that);
    int get_size () const;
	T get_determinant()  const;
    SquareMatrix get_inverse() const;
    SquareMatrix get_trans () const;
    SquareMatrix& transpose ();
    SquareMatrix<T>& invert();
    T get_trace() const;
    SquareMatrix operator-(const SquareMatrix& ) const;
    SquareMatrix operator+(const SquareMatrix& ) const;
    SquareMatrix operator*(const SquareMatrix& ) const;
    SquareMatrix& operator=(const SquareMatrix<T>& );
    SquareMatrix operator*(const T k) const;
    SquareMatrix& operator+=(const SquareMatrix& );
    SquareMatrix& operator-=(const SquareMatrix& );
    SquareMatrix& operator*=(const SquareMatrix& );
    SquareMatrix& operator*=(const T k);
    template <typename C>
    friend Matrix<C> operator*(const C k, const Matrix<C>& a);
};


template <typename T>
void SquareMatrix<T>::assign(const Matrix<T>& that) {
    Matrix<T>::assign(that);
}
template <typename T>
SquareMatrix<T>::SquareMatrix(int size){
    Matrix<T>::set_mem(size, size);
    for (int i = 0; i < Matrix<T>::get_rows_numb(); i++) {
        for (int j = 0; j < Matrix<T>::get_columns_numb(); j++) {
            Matrix<T>::set (i, j, 0);
        }
    }
}
template <typename T>
SquareMatrix<T>::SquareMatrix(){}

template <typename T>
SquareMatrix<T>::SquareMatrix(const Matrix<T>& that){
    assign(that);
}

template <typename T>
T SquareMatrix<T>::get_determinant() const {
    SquareMatrix<T> temp(*this);
    bool f = true;
    T det = T(1);
    for (int col = 0; col < temp.get_size() && f == true; col++) {
        int j = col;
        while (j < temp.get_size() && temp.Matrix<T>::get(j,col) == T(0)) {
            j++;
        }
        if (j == temp.get_size()) {
            f = false;
        }
        else {
            if (j != col) {
                temp.Matrix<T>::swap_rows(col, j);
            }
            for (int k = col + 1; k < temp.get_size(); k++)
            {
                T coef = temp.Matrix<T>::get(k, col) / temp.Matrix<T>::get(col, col);
                temp.Matrix<T>::row_add(k, col, -coef);
            }
        }
        det *= temp.Matrix<T>::get(col, col);
    }
    if (!f) {
        det = T(0);
    }
    return det; 
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::get_minor(const int row,const int col) const {
    SquareMatrix res;
    res.Matrix<T>::set_mem(get_size() - 1, get_size() - 1);
    for (int i = 0; i < res.get_size(); i++) {
        for (int j = 0; j < res.get_size(); j++) {
            int di = 0, dj = 0;
            if (i >= row)
                di = 1;
            if (j >= col)
                dj = 1;
            res.Matrix<T>::set(i, j, this->Matrix<T>::get(i + di, j + dj));
        }
    }
    return res;
}
template <typename T>
SquareMatrix<T> SquareMatrix<T>::get_trans () const {
    SquareMatrix res(get_size());
    for (int i = 0; i < get_size(); i++) {
        for (int j = 0; j < get_size(); j++) {
            res.Matrix<T>::set(i, j, this->Matrix<T>::get(j, i));
        }
    }
    return res;
}
template <typename T>
SquareMatrix<T>& SquareMatrix<T>::transpose() {
    *this = get_trans();
    return *this;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::get_inverse() const {
    if (this->get_determinant() == T(0)) {
        throw MatrixIsDegenerateError();
    }
    SquareMatrix res;
    res.Matrix<T>::set_mem(get_size(), get_size());
    for (int i = 0; i < res.get_size(); i++) {
        for (int j = 0; j < res.get_size(); j++) {
            T a = 0;
            if ((i + j) % 2 == 0) {
                a =  get_minor(i, j).get_determinant() ;
            } else {
                a = -get_minor(i, j).get_determinant() ;
            }
            res.Matrix<T>::set(i, j, a);
        }
    }
    res = res.transpose();
    res *= (1 / get_determinant());
    return res;
}
template <typename T>
SquareMatrix<T>& SquareMatrix<T>::invert() {
    *this = get_inverse();
    return *this;
}
template <typename T>
T SquareMatrix<T>::get_trace() const {
    T sum = 0;
    for (int i = 0; i < this->get_size(); i++) {
        sum += this->Matrix<T>::get(i, i);
    }
    return sum;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator+(const SquareMatrix<T>& that) const {
    Matrix<T> newthis(*this);
    Matrix<T> newthat(that);
    SquareMatrix res(newthis + newthat);
    return res;
}
template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator-(const SquareMatrix<T>& that) const {
    return SquareMatrix(*this + that * (-1));
}
template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(const SquareMatrix<T>& that) const {
    Matrix<T> newthis(*this);
    Matrix<T> newthat(that);
    SquareMatrix res(newthis * newthat);
    return res;
}
template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator=(const SquareMatrix<T>& that) {
    if (this != &that) {
        Matrix<T>::delete_mem();
        assign(that);
    }
    return *this;
}
template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator+=(const SquareMatrix& that) {
    *this = *this + that;
    return *this;
}
template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator-=(const SquareMatrix& that) {
    *this = *this - that;
    return *this;
}
template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(const SquareMatrix& that) {
    *this = *this * that;
    return *this;
}
template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(const T k) {
    *this = *this * k;
    return *this;
}
template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(const T k) const {
    Matrix<T> newthis(*this);
    SquareMatrix res(newthis * k);
    return res;
}
template <typename C>
SquareMatrix<C> operator*(const C k, const SquareMatrix<C>& a) {
    return a * k;
}
template <typename T>
int SquareMatrix<T>::get_size() const {
    return this->Matrix<T>::get_rows_numb();
}


//Main

int main() {
    int m, n, p, q;
    cin >> m >> n >> p >> q;
    
    Matrix<int> A(m, n), B(p, q);
    cin >> A >> B;
    
    A = A;
    try {
        cout << A + B * 2 - m * A << endl;
        cout << (A -= B += A *= 2) << endl;
        cout << (((A -= B) += A) *= 2) << endl;
    } catch (const MatrixWrongSizeError&) {
        cout << "A and B are of different size." << endl;
    }
    B = A;
    cout << B << endl;
    
    Rational r;
    cin >> r;
    Matrix<Rational> C(m, n), D(p, q);
    cin >> C >> D;
    try {
        cout << C * D << endl;
        cout << (C *= D) << endl;
        cout << C << endl;
    } catch (const MatrixWrongSizeError&) {
        cout << "C and D have not appropriate sizes for multiplication." << endl;
    }
    cout << C.get_trans() * (r * C) << endl;
    cout << C.transpose() << endl;
    cout << C << endl;
    
    SquareMatrix<Rational> S(m);
    cin >> S;
    SquareMatrix<Rational> P(S);
    const SquareMatrix<Rational>& rS = S;
    cout << rS.get_size() << ' ' << rS.get_determinant() << ' ' << rS.get_trace() << endl;
    cout << (S = S) * (S + rS) << endl;
    cout << (S *= S) << endl;
    C.transpose();
    cout << rS * C << endl;
    cout << S << endl;
    S = P;
    cout << (Rational(1, 2) * S).get_determinant() << endl;
    try {
        cout << rS(0, 0) << endl;
        (S(0, 0) *= 2) /= 2;
        cout << rS(0, m) << endl;
    } catch (const MatrixIndexError&) {
        cout << "Index out of range." << endl;
    }
    cout << rS.get_trans() << endl;
    try {
        cout << rS.get_inverse() << endl;
        cout << S.invert().get_trans().get_determinant() << endl;
        cout << S << endl;
    } catch (const MatrixIsDegenerateError&) {
        cout << "Cannot inverse S." << endl;
    } 
    return 0;
}