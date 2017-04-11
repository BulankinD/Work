#include<iostream>
#include<cstring>
#include<cstdlib>
#include<stdexcept>
using namespace std;

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

int main(int argc, char** argv) {
    int a;
    cin >> a;

    int p, q;
    cin >> p >> q;
    const Rational rc(p, q); // q != 0 is guaranteed by author of tests
    cout << rc.getNumerator() << ' ' << rc.getDenominator() << endl;

    Rational r1, r2;
    cin >> r1 >> r2;

    cout << r1 << endl;
    cout << r2 << endl;

    try {
        cout << 1/r1 << endl;
    } catch (const RationalDivisionByZero& ex) {
        cout << "Cannot get reciprocal of r1." << endl;
    }

    try {
        cout << rc/r2 << endl;
    } catch (const RationalDivisionByZero& ex) {
        cout << "Cannot divide by r2." << endl;
    }

    cout << (r1 < r2) << endl;
    cout << (r1 <= r2) << endl;
    cout << (r1 > r2) << endl;
    cout << (r1 >= r2) << endl;
    cout << (r1 == r2) << endl;
    cout << (r1 != r2) << endl;

    cout << (r1 < a) << endl;
    cout << (r1 <= a) << endl;
    cout << (r1 > a) << endl;
    cout << (r1 >= a) << endl;
    cout << (r1 == a) << endl;
    cout << (r1 != a) << endl;

    cout << (a < r2) << endl;
    cout << (a <= r2) << endl;
    cout << (a > r2) << endl;
    cout << (a >= r2) << endl;
    cout << (a == r2) << endl;
    cout << (a != r2) << endl;

    cout << rc + a << endl
         << a + rc << endl
         << -rc * r1 << endl
         << (+r1 - r2 * rc) * a << endl;

    cout << ++r1 << endl;
    cout << r1 << endl;
    cout << r1++ << endl;
    cout << r1 << endl;
    cout << --r1 << endl;
    cout << r1 << endl;
    cout << r1-- << endl;
    cout << r1 << endl;
    cout << ++++r1 << endl;
    cout << r1 << endl;

    cout << ((((r1 += r2) /= Rational(-5,3)) -= rc) *= a) << endl;
    cout << (r1 += r2 /= 3) << endl;
    cout << r1 << endl;
    cout << r2 << endl;
    return 0;
}