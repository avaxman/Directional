// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>

//This header file implements a "home-made" Big Integer type, which is only needed in case GMP is not installed. Note: this is slow-ish.

class BigInteger {
private:
    static const long long BASE = 1e9;
    static const int DIGITS = 9;
    static const int CONVERTIBLE_SIZE = 2;   //maximum amount of digits in an integer that can be safely back-converted to long longs
    
    std::vector<long long> digits;
    bool negative;
    
    void trim() {
        int whileTest=0;
        while (!digits.empty() && digits.back() == 0) {
            digits.pop_back();
            whileTest++;
            assert(whileTest<10000 && "trim(): while running too long! ");
        }
        if (digits.empty()) {
            digits.push_back(0);
            negative = false;
        }
    }
    
    
    
public:
    BigInteger() : negative(false) {}
    
    BigInteger(long long value) {
        if (value < 0) {
            negative = true;
            value = -value;
        } else {
            negative = false;
        }
        
        digits.clear();
        if (value==0){
            digits.push_back(0);
            return;
        }
        while (value > 0) {
            digits.push_back(value % BASE);
            value /= BASE;
        }
    }
    
    long long convert() const{
        assert(digits.size()<=CONVERTIBLE_SIZE) && "integer is not of a convertible size!");
        long long result=0.0;
        for (int i=digits.size()-1; i>=0 ;i--)
            result = result * BASE + digits[i];
        
        return (negative ? -result : result);
    }
    
    BigInteger abs() const{
        BigInteger absNum = *this;
        absNum.negative=false;
        return absNum;
    }
    
    BigInteger operator+(const BigInteger &other) const {
        if (other==0)
            return *this;
        if (*this==0)
            return other;
        if (negative != other.negative) {
            return *this - (-other);
        }
        
        BigInteger result;  //TODO: reserve
        result.negative = negative;
        
        long long carry = 0;
        for (size_t i = 0; (i < std::max(digits.size(), other.digits.size())) || carry; i++) {
            long long sum = carry;
            if (i < digits.size()) sum += digits[i];
            if (i < other.digits.size()) sum += other.digits[i];
            
            result.digits.push_back(sum % BASE);
            carry = sum / BASE;
        }
        
        result.trim();
        return result;
    }
    
    BigInteger operator-() const {
        BigInteger result = *this;
        if (result!=0)
            result.negative = !negative;
        
        return result;
    }
    
    BigInteger operator-(const BigInteger &other) const {
        if (other==0)
            return *this;
        if (*this==0)
            return -other;
        if (negative != other.negative) {
            return *this + (-other);
        }
        
        if (this->abs() < other.abs()) {
            return -(other - *this);
        }
        
        BigInteger result;   //TODO: reserve
        result.negative = negative;
        
        long long borrow = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            long long diff = digits[i] - borrow;
            if (i < other.digits.size()) diff -= other.digits[i];
            
            if (diff < 0) {
                diff += BASE;
                borrow = 1;
            } else {
                borrow = 0;
            }
            
            result.digits.push_back(diff);
        }
        
        result.trim();
        return result;
    }
    
    BigInteger operator*(const BigInteger &other) const {
        BigInteger result;
        result.digits.resize(digits.size() + other.digits.size());
        result.negative = (negative != other.negative);
        
        for (size_t i = 0; i < digits.size(); ++i) {
            long long carry = 0;
            for (int j = 0; (j < other.digits.size()) || carry; ++j) {
                long long prod = result.digits[i + j] + digits[i] * (j < other.digits.size() ? other.digits[j] : 0) + carry;
                result.digits[i + j] = prod % BASE;
                carry = prod / BASE;
            }
        }
        
        result.trim();
        return result;
    }
    
    //Specialized version for multiplying directly with a single digit (effective for division)s
    BigInteger operator*(const long long &other) const {
        BigInteger result;
        result.digits.resize(digits.size() + 1);
        result.negative = (negative != other<0);
        
        long long carry = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            long long prod = digits[i] * other + carry;
            result.digits[i] = prod % BASE;
            carry = prod / BASE;
        }
        result.digits[digits.size()]=carry;
        result.trim();
        return result;
    }
    
    //A version that expects the result to be a single digit (like as a part of long division).
    inline long long single_digit_division(const BigInteger& other, BigInteger& mod) const{
        
        if (*this==0){
            mod = 0;
            return 0;
        }
        
        if (other==1){
            mod = 0;
            return this->convert();
        }
        
        if (this->abs() < other.abs()){
            mod = *this;
            return 0;
        }
        
        
        if ((this->digits.size()<=CONVERTIBLE_SIZE)&&(other.digits.size()<=CONVERTIBLE_SIZE)){
            long long convertThis = this->convert();
            long long convertOther = other.convert();
            long long quotient = convertThis/convertOther;
            mod = BigInteger(convertThis - convertOther*quotient);   //TODO: is this corret and not overflowing?
            return quotient;
        }
        
        //Evaluating the result by dividing digit by digit and refining the result
        BigInteger dividend = this->abs();
        BigInteger divisor = other.abs();
        
        BigInteger currDividend, currDivisor;
        currDividend.digits.resize(1+dividend.digits.size()-divisor.digits.size());
        currDivisor.digits.resize(1);
        currDivisor.digits[0]=divisor.digits[divisor.digits.size()-1];
        for (int i=0;i<currDividend.digits.size();i++)
            currDividend.digits[currDividend.digits.size()-i-1] = dividend.digits[dividend.digits.size()-1-i];
        
        long long quotient = currDividend.convert()/currDivisor.convert();
        long long left = (currDividend).convert()/(currDivisor+1).convert();
        long long right = (currDividend+1).convert()/(currDivisor).convert();

        int whileTest=0;
        BigInteger diff; diff.digits.reserve(dividend.digits.size());
        while (left <= right) {
            long long mid = (left + right) / 2;
            diff = dividend - divisor * mid;
            if (diff >= 0) {
                quotient = mid;
                mod = diff;
                if (diff==0)
                    break; //it's found
                left = mid + 1;
            } else {
                right = mid - 1;
            }
            whileTest++;
            assert("operator/: while running too long! " && whileTest<10000);
        }

        return(this->negative != other.negative ? -quotient : quotient);
    }
    
    
    
    inline BigInteger operator/(const BigInteger &other) const {
        
        //rule out simple cases
        if (*this==0)
            return 0;
        
        if (other==1)
            return *this;
        
        if (this->abs() < other.abs())
            return 0;
        
        if ((this->digits.size()<=CONVERTIBLE_SIZE)&&(other.digits.size()<=CONVERTIBLE_SIZE)){
            long long convertThis = this->convert();
            long long convertOther = other.convert();
            return convertThis/convertOther;
        }
        
        
        //TODO: leading zeros
        
        BigInteger dividend = this->abs();
        BigInteger divisor = other.abs();
        
        BigInteger quotient;
        quotient.digits.resize(dividend.digits.size());
        BigInteger current;
        BigInteger mod;
        for (int i = dividend.digits.size()-1; i>=0 ; i--) {
            current.digits.insert(current.digits.begin(), dividend.digits[i]);
            quotient.digits[i] = current.single_digit_division(divisor, mod);  //updates current as the modulo
            current = mod;
            
        }
        
        quotient.negative = (negative != other.negative);
        quotient.trim();

        return quotient;
    }
    
    BigInteger operator%(const BigInteger &other) const {
        return *this - (*this / other) * other;
    }
    
    BigInteger &operator/=(const BigInteger &other) {
        *this = *this / other;
        return *this;
    }
    
    bool operator<(const BigInteger &other) const {
        if (negative != other.negative) {
            return negative;
        }
        
        if (digits.size() != other.digits.size()) {
            return (digits.size() < other.digits.size()) ^ negative;
        }
        
        for (size_t i = digits.size(); i-- > 0;) {
            if (digits[i] != other.digits[i]) {
                return (digits[i] < other.digits[i]) ^ negative;
            }
        }
        
        return false;
    }
    
    bool operator>(const BigInteger &other) const {
        return other < *this;
    }
    
    bool operator<=(const BigInteger &other) const {
        return *this < other || *this == other;
    }
    
    bool operator>=(const BigInteger &other) const {
        return !(*this < other);
    }
    
    bool operator==(const BigInteger &other) const {
        return digits == other.digits && negative == other.negative;
    }
    
    bool operator!=(const BigInteger &other) const {
        return !(*this == other);
    }
    
    BigInteger operator++() {
        *this = *this + BigInteger(1);
        return *this;
    }
    
    std::string to_string() const {
        std::ostringstream oss;
        if (negative) oss << '-';
        //if (digits.empty()) {
        //    oss << '0';
        //} else {
        oss << digits.back();
        for (size_t i = digits.size() - 1; i-- > 0;) {
            oss << std::setw(DIGITS) << std::setfill('0') << digits[i];
        }
        //}
        return oss.str();
    }
    
    friend std::ostream &operator<<(std::ostream &out, const BigInteger &value) {
        out << value.to_string();
        return out;
    }
};

BigInteger gcd(BigInteger a, BigInteger b) {
    int whileTest=0;
    if (b > a) std::swap(a,b);
    if ((b==1)||(a==1))
        return 1;
    while (b != BigInteger(0)) {
        BigInteger temp = b;
        b = a % b;
        a = temp;
        whileTest++;
        assert(whileTest<10000 && "gcd(): while running too long! ");
    }
    return a;
}
