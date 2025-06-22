#ifndef DIRECTIONAL_ENUMBER_GMP_HEADER_FILE
#define DIRECTIONAL_ENUMBER_GMP_HEADER_FILE

#include <gmpxx.h>

class EInt {
public:
    mpz_class value;

    EInt() : value(0) {}
    EInt(long long v) : value(static_cast<signed long>(v)) {}
    EInt(const mpz_class& v) : value(v) {}
    EInt(const EInt& other) = default;

    EInt& operator=(const EInt& other) = default;

    EInt abs() const { return EInt(::abs(value)); }

    long long convert() const {
        return value.get_si();  // only if fits in long long
    }

    std::string to_string() const {
        return value.get_str();
    }

    EInt operator+(const EInt& other) const { return EInt(value + other.value); }
    EInt operator-(const EInt& other) const { return EInt(value - other.value); }
    EInt operator-() const { return EInt(-value); }
    EInt operator*(const EInt& other) const { return EInt(value * other.value); }
    EInt operator/(const EInt& other) const { return EInt(value / other.value); }
    EInt operator%(const EInt& other) const { return EInt(value % other.value); }

    EInt& operator+=(const EInt& other) { value += other.value; return *this; }
    EInt& operator-=(const EInt& other) { value -= other.value; return *this; }
    EInt& operator*=(const EInt& other) { value *= other.value; return *this; }
    EInt& operator/=(const EInt& other) { value /= other.value; return *this; }

    bool operator==(const EInt& other) const { return value == other.value; }
    bool operator!=(const EInt& other) const { return value != other.value; }
    bool operator<(const EInt& other) const { return value < other.value; }
    bool operator<=(const EInt& other) const { return value <= other.value; }
    bool operator>(const EInt& other) const { return value > other.value; }
    bool operator>=(const EInt& other) const { return value >= other.value; }

    friend std::ostream& operator<<(std::ostream& os, const EInt& ei) {
        os << ei.value;
        return os;
    }

    const mpz_class& raw() const { return value; }
};

inline EInt gcd(EInt a, EInt b) {
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), a.value.get_mpz_t(), b.value.get_mpz_t());
    return EInt(g);
}

class ENumber {
public:
    mpq_class value;

    ENumber() : value(0) {}
    ENumber(double x, double tol = 1e-9) : value(x) {}

    ENumber(const EInt& _num, const EInt& _den, bool simplify = true) {
        value = mpq_class(_num.value, _den.value);
        if (simplify) value.canonicalize();
    }

    ENumber(const EInt& _num) : value(mpq_class(_num.value)) {}
    ENumber(const mpq_class& val) : value(val) {}

    ENumber& operator=(const ENumber& other) = default;

    ENumber operator+(const ENumber& b) const { return ENumber(value + b.value); }
    ENumber operator-(const ENumber& b) const { return ENumber(value - b.value); }
    ENumber operator-() const { return ENumber(-value); }
    ENumber operator*(const ENumber& b) const { return ENumber(value * b.value); }
    ENumber operator/(const ENumber& b) const { return ENumber(value / b.value); }

    ENumber& operator+=(const ENumber& b) { value += b.value; return *this; }
    ENumber& operator/=(const ENumber& b) { value /= b.value; return *this; }

    bool operator==(const ENumber& b) const { return value == b.value; }
    bool operator!=(const ENumber& b) const { return value != b.value; }
    bool operator<(const ENumber& b) const { return value < b.value; }
    bool operator<=(const ENumber& b) const { return value <= b.value; }
    bool operator>(const ENumber& b) const { return value > b.value; }
    bool operator>=(const ENumber& b) const { return value >= b.value; }

    ENumber abs() const { return ENumber(::abs(value)); }

    long double to_double(int maxDigits = 12) const {
        return value.get_d();
    }

    EInt num() const { return EInt(value.get_num()); }
    EInt den() const { return EInt(value.get_den()); }

    const mpq_class& raw() const { return value; }
};

#endif
