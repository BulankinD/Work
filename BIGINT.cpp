#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <functional>

int step(int, int);

class BigInt {
private:
	int sign_;
	std::vector<unsigned int> data_;

    static const unsigned int BASE_LG = 9;
	static const unsigned long long int BASE = 1000000000;

	void remove_lead_zeros();
    

public:
    BigInt(int value);
	BigInt(long int value);
	BigInt(unsigned int value);
	BigInt(unsigned long long int value);
	BigInt(std::string value);
	BigInt(const BigInt&);
	BigInt();

    BigInt operator-() const;
    BigInt& operator++();
    BigInt& operator--();
    const BigInt operator++(int);
    const BigInt operator--(int);

    BigInt& operator=(const BigInt& other);

    BigInt& operator-=(const BigInt& other);
    BigInt& operator+=(const BigInt& other);
    BigInt& operator*=(const BigInt& other);
    BigInt& operator/=(const BigInt& other);
    BigInt& operator%=(const BigInt& other);

    std::string toString() const;
	BigInt& div2();
	BigInt& shl(size_t cnt);

    BigInt operator+(const BigInt& other) const;
	BigInt operator+(const int& other) const;
    BigInt operator-(const BigInt& other) const;
    BigInt operator*(const BigInt& other) const;
	BigInt operator*(const int& other) const;
    BigInt operator/(const BigInt& other) const;
	BigInt operator/(const int other) const;
    BigInt operator%(const BigInt& other) const;
	BigInt sqrt() const;

    bool operator<(const BigInt& other) const;
    bool operator>(const BigInt& other) const;
    bool operator<=(const BigInt& other) const;
    bool operator>=(const BigInt& other) const;
    bool operator==(const BigInt& other) const;
    bool operator!=(const BigInt& other) const;

    friend std::istream& operator>>(std::istream& is, BigInt& big);
};

std::ostream& operator<<(std::ostream& os, const BigInt& big);
std::string to_string(unsigned long int that);



BigInt::BigInt() {
	this->sign_ = 0;
}

BigInt nod(const BigInt& a, const BigInt& b) {
	if (a == 0 || b == 0)
		return a + b;
	if (a > b) 
		return nod(b, a % b);
	else 
		return nod(a, b % a);
}

BigInt::BigInt(const BigInt& value) {
	this->sign_ = value.sign_;
	for (size_t i = 0; i < value.data_.size(); i++)
		data_.push_back(value.data_[i]); 
}

BigInt::BigInt(int value) {
    sign_ = (0 < value) - (value < 0);
    value *= sign_;
    while (value != 0) {
        data_.push_back(value % BASE);
        value /= BASE;	
    }
}

BigInt::BigInt(unsigned int value) {
    sign_ = 1;
    while (value != 0) {
        data_.push_back(value % BASE);
        value /= BASE;
    }
}

BigInt::BigInt(long int value) {
    sign_ = (0 < value) - (value < 0);
    value *= sign_;
    while (value != 0) {
        data_.push_back(value % BASE);
        value /= BASE;
    }
}

BigInt::BigInt(unsigned long long int value) {
    sign_ = 1;
    while (value != 0) {
        data_.push_back(value % BASE);
        value /= BASE;
    }
}

BigInt::BigInt(std::string value) {
	if (value[0] == '-') {
		sign_ = -1;
		int len = value.length();
		int n = (len - 1) / BASE_LG;
		int l = (len - 1) % BASE_LG;
		for (int i = 0; i < n; i++) {
			int sum = 0;
			for (int j = 0; j < BASE_LG; j++) {
				sum += (int(value[len - 1 - i * BASE_LG - j]) - 48) * step(10, j);
			}
			data_.push_back(sum);
		}
		int sum = 0;
		for (int i = l - 1; i > 0; i--)
			sum += (int(value[i]) - 48) * step(10, l - 1 - i);
		if (sum != 0)
			data_.push_back(sum);
	} else if (value[0] == '0') {
		sign_ = 0;
	} else { 
		sign_ = 1;
		int len = value.length();
		int n = len / BASE_LG;
		int l = len % BASE_LG;
		for (int i = 0; i < n; i++) {
			int sum = 0;
			for (int j = 0; j < BASE_LG; j++) {
				sum += (value[len - 1 - i * BASE_LG - j] - 48) * step(10, j);
			}
			data_.push_back(sum);
		}
		int sum = 0;
		for (int i = l - 1; i >= 0; i--) {
			sum += (int(value[i]) - 48) * step(10, l - 1 - i);
		}
		if (sum != 0)
			data_.push_back(sum);
    }
}

void BigInt::remove_lead_zeros() {
    while (data_.size() && !data_.back())
        data_.pop_back();
    sign_ = data_.size() ? sign_ : 0;
}

BigInt& BigInt::div2() {
    for (size_t i = data_.size() - 1; i != static_cast<size_t>(-1); --i) {
        if (i != 0) 
            data_[i - 1] += (data_[i] & 1) * BASE;
        data_[i] >>= 1;
    }
    remove_lead_zeros();
    return *this;
}

BigInt BigInt::operator-() const {
    BigInt res = *this;
    res.sign_ = -res.sign_;
    return res;
}

BigInt& BigInt::operator++() {
    return *this += 1;
}

BigInt& BigInt::operator--() {
    return *this -= 1;
}

const BigInt BigInt::operator++(int) {
    BigInt res = *this;
    ++(*this);
    return res;
}

const BigInt BigInt::operator--(int) {
    BigInt res = *this;
    --(*this);
    return res;
}

BigInt& BigInt::operator=(const BigInt& other) {
	this->sign_ = other.sign_;
	this->data_.clear();
	for (size_t i = 0; i < other.data_.size(); i++)
		this->data_.push_back(other.data_[i]);
	return *this;
}

BigInt BigInt::operator+(const BigInt& other) const {
    BigInt res(*this);
    return res += other;
}

BigInt BigInt::operator+(const int& other) const {
    BigInt a(other);
    return *this + a;
}

BigInt BigInt::operator-(const BigInt& other) const {
    BigInt res(*this);
    return res -= other;
}

BigInt BigInt::operator*(const BigInt& other) const {
    BigInt res(*this);
    return res *= other;
}

BigInt BigInt::operator*(const int& other) const {
    BigInt a(other);
    return *this * a;
}

BigInt BigInt::operator/(const BigInt& other) const {
    BigInt res(*this);
    return res /= other;
}

BigInt BigInt::operator/(const int other) const {
    BigInt a(other);
    return *this / a;
}

BigInt BigInt::operator%(const BigInt& other) const {
    BigInt res(*this);
    return res %= other;
}

std::string BigInt::toString() const {
    std::string res = sign_ < 0 ? "-" : "";
    res += sign_ ? "" : "0";
	if (data_.size() != 0)
		res += to_string(data_[data_.size() - 1]);
    for (int i = data_.size() - 2; i >= 0; i--) {
        std::string digit = to_string(data_[i]);
		int k = BASE_LG - digit.length();
		std::string tmp = "";
		for (int j = 0; j < k; j++)
			tmp += "0";
        res += (tmp + digit);
	}
    return res;
}

BigInt& BigInt::operator+=(const BigInt& other) {
    if ((sign_ < 0) ^ (other.sign_ <= 0))
        return *this -= -other;

    unsigned int transfer = 0;
    for (size_t i = 0; i < other.data_.size() || transfer; ++i) {
        if (i == data_.size())
            data_.push_back(0);

        data_[i] += (other.data_.size() > i ? other.data_[i] : 0) + transfer;
        transfer = data_[i] >= BASE;
        data_[i] -= transfer * BASE;
    }
    sign_ = sign_ ? sign_ : other.sign_;
    remove_lead_zeros();
    return *this;
}

BigInt& BigInt::operator-=(const BigInt& other) {
    if ((sign_ < 0) ^ (other.sign_ < 0))
        return *this += -other;
    unsigned int loan = 0;
    for (size_t i = 0; i < other.data_.size() || (i < data_.size() && loan); ++i) {
        if (i == data_.size())
            data_.push_back(0);

        unsigned int subtrahend = (other.data_.size() > i ? other.data_[i] : 0) + loan;
        loan = data_[i] < subtrahend;
        data_[i] += loan * BASE;
        data_[i] -= subtrahend;
    }
    if (loan) {
        sign_ = -other.sign_;
        unsigned int transfer = 0;
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] = BASE - (i != 0) - data_[i] + transfer;
            transfer = data_[i] >= BASE;
            data_[i] -= transfer * BASE;
        }
    }
    remove_lead_zeros();
    return *this;
}

BigInt& BigInt::operator*=(const BigInt& other) {
    sign_ *= other.sign_;
    BigInt bufferA(*this);
    BigInt bufferB(other);
    data_.assign(data_.size() + bufferB.data_.size(), 0);
    for (size_t i = 0; i < bufferA.data_.size(); ++i) {
        unsigned long long int transfer = 0;
        for (size_t j = 0; j < bufferB.data_.size() || transfer; ++j) {
            unsigned long long int digit = data_[i + j] + transfer + 
                static_cast<unsigned long long int>(bufferA.data_[i]) * (j < bufferB.data_.size() ? bufferB.data_[j] : 0);

            data_[i + j] = digit % BASE;
            transfer = digit / BASE;
        }
    } 
    remove_lead_zeros();
    return *this;
}

BigInt& BigInt::shl(size_t cnt) {
    size_t old_size = data_.size();
    data_.resize(data_.size() + cnt);
    for (size_t i = old_size - 1; i != static_cast<size_t>(-1); --i)
        data_[i + cnt] = data_[i];
    for (size_t i = 0; i < cnt; ++i)
        data_[i] = 0;

    return *this;
}

BigInt& BigInt::operator/=(const BigInt& other) {
    if (sign_ * other.sign_ == -1) {
        BigInt tmp(*this * sign_);
        tmp /= (other * other.sign_);
        return *this = (tmp *= -1);
    }

    if (sign_ == -1 && other.sign_ == -1) {
        BigInt tmp(*this * sign_);
        return *this = (tmp /= (other * other.sign_));
    }

    if (*this < other)
        return *this = 0;

    if (*this == other)
        return *this = 1;

    bool high_comp = data_.back() < other.data_.back() + 1;
    int cnt = data_.size() - (other.data_.size() + high_comp);
    if (cnt < 0)
        return *this = 1;
    cnt = cnt > 0 ? cnt : 0;
    unsigned long long int dividend = data_.back();
    if (high_comp) {
        dividend *= BASE;
        dividend += data_[data_.size() - 2];
    }

    unsigned int divider = other.data_.back() + 1;
    unsigned long long int quotient = dividend / divider;
    BigInt res(quotient);
    BigInt tmp = res * other;
    tmp.shl(cnt);
    return *this = res.shl(cnt) + ((*this - tmp) / other);
}

BigInt& BigInt::operator%=(const BigInt& other) {
    *this -= (*this / other) * other;
    return *this;
}

bool BigInt::operator<(const BigInt& other) const {
    if (sign_ != other.sign_) 
        return sign_ < other.sign_;      
    if (data_.size() != other.data_.size())
        return (data_.size() < other.data_.size()) ^ (sign_ < 0);
    for (size_t i = data_.size() - 1; i != static_cast<size_t>(-1); --i)
        if (data_[i] != other.data_[i])
            return (data_[i] < other.data_[i]) ^ (sign_ < 0);

    return 0;
}

bool BigInt::operator>(const BigInt& other) const {
    return other < *this;
}

bool BigInt::operator<=(const BigInt& other) const {
    return !(other < *this);
}

bool BigInt::operator>=(const BigInt& other) const {
    return !(*this < other);
}

bool BigInt::operator==(const BigInt& other) const {
    return !(*this < other) && !(other < *this);
}

bool BigInt::operator!=(const BigInt& other) const {
    return (*this < other) || (other < *this);
}

BigInt BigInt::sqrt() const {
    BigInt r(*this);
    BigInt l(1);
    BigInt mid((r + l).div2());
    BigInt ans(mid * mid);
    bool f = true;
    while (f) {
        if (ans == *this || (*this > ans && (mid + 1) * (mid + 1) > *this)) {
            f = false;
        }
        else {
            if (*this > ans) {
                l = mid;
            } else {
                r = mid;
            }
            mid = (r + l).div2();
            ans = mid * mid;
        }
    }
    return mid;
}

std::ostream& operator<<(std::ostream& os, const BigInt& big) {
    os << big.toString();
    return os;
}

std::istream& operator>>(std::istream& is, BigInt& big) {
    std::string sdata;
    is >> sdata;
	BigInt ans(sdata);
	big = ans;
    return is;
}

std::string to_string(unsigned long int that) {
	unsigned int a = that;
	std::string str = "";
	while (a > 0) {
		str += char(a % 10 + 48);
		a /= 10;
	}
	std::string ans = "";
	for(int i = str.length() - 1; i >= 0; i--)
		ans += str[i];
	return ans;
}

int main() {	
}

int step(int a, int b) {
	int ans = 1;
	for (int i = 0; i < b; i++)
		ans *= a;
	return ans;
}