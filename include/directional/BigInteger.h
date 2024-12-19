#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>

class BigInteger {
private:
    static const long long BASE = 1e9;
    static const int DIGITS = 9;

    std::vector<long long> digits;
    bool negative;

    void trim() {
        int whileTest=0;
        while (!digits.empty() && digits.back() == 0) {
            digits.pop_back();
            whileTest++;
            assert("trim(): while running too long! " && whileTest<10000);
        }
        if (digits.empty()) {
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
        while (value > 0) {
            digits.push_back(value % BASE);
            value /= BASE;
        }
    }
  
  BigInteger abs() const{
    BigInteger absNum = *this;
    absNum.negative=false;
    return absNum;
  }

    /*BigInteger(const std::string &str) {
        negative = false;
        digits.clear();

        int start = 0;
        if (str[0] == '-') {
            negative = true;
            start = 1;
        }

        for (int i = str.length(); i > start; i -= DIGITS) {
            int end = i;
            int begin = std::max(start, i - DIGITS);
            digits.push_back(std::stoll(str.substr(begin, end - begin)));
        }

        trim();
    }*/

    BigInteger operator+(const BigInteger &other) const {
        if (other==0)
            return *this;
        if (*this==0)
            return other;
        if (negative != other.negative) {
            return *this - (-other);
        }

        BigInteger result;
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
        if (!result.digits.empty()) {
            result.negative = !negative;
        }
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

        BigInteger result;
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
      /*std::cout<<"this: "<<this->to_string()<<std::endl;
      std::cout<<"other: "<<other.to_string()<<std::endl;
      std::cout<<"product: "<<result.to_string()<<std::endl;*/
      return result;
    }

    /*BigInteger operator/(const BigInteger &other) const {
        assert("Biginteger division: division by zero!" && other != BigInteger(0));
        BigInteger dividend = *this;
        BigInteger divisor = other;
        dividend.negative = divisor.negative = false;

        BigInteger quotient;
        quotient.digits.resize(dividend.digits.size());
        BigInteger currDividend;  //the current part of the dividend we are trying to fit in
        for (int i = dividend.digits.size()-1 ; i>=0 ; i--){
            currDividend = currDividend * BigInteger(BASE) + BigInteger(dividend.digits[i]);
            if (divisor < currDividend){
                quotient = quotient*BigInteger(BASE);
                continue;
            }

            //finding out how many times the divisor fits in the currDividend
            //doing binary search
            long long currQuotient = BASE-1;
            long long currLeap = BASE/2;
            for (int currLeap=BASE/2 ; currLeap >= 1 ; currLeap/=2){
                BigInteger diff = BigInteger(currQuotient)*divisor - BigInteger(currDividend);
                if (diff>0)
                    currQuotient-=currLeap;
                if (diff<0) 
                    currQuotient+=currLeap; 
                if (diff==0)
                    break;  //exactly multiplies
            }

            quotient = quotient*BigInteger(BASE)+currQuotient;
            currDividend = currDividend - BigInteger(currQuotient)*divisor;
        }

    }*/
    
    
    
    BigInteger operator/(const BigInteger &other) const {
        

        BigInteger dividend = *this;
        BigInteger divisor = other;
        dividend.negative = divisor.negative = false;

        if (divisor==1)
          int kaka=8;
        BigInteger quotient;
        quotient.digits.resize(dividend.digits.size());
        //std::cout<<"Dividing "<<dividend.to_string()<<" by "<<divisor.to_string()<<std::endl;
        BigInteger current;
        for (int i = dividend.digits.size()-1; i>=0 ; i--) {
            current = current * BigInteger(BASE) + BigInteger(dividend.digits[i]);
            //std::cout<<"current: "<<current.to_string()<<std::endl;
            long long count = 0;
            long long left = 0, right = BASE - 1;
            int whileTest=0;
            while (left <= right) {
                //std::cout<<"left: "<<left<<" right:"<<right<<std::endl;
                long long mid = (left + right) / 2;
                //std::cout<<"mid: "<<mid<<std::endl;
                //std::cout<<"divisor * BigInteger(mid): "<<(divisor * BigInteger(mid)).to_string()<<std::endl;
                BigInteger diff = current - divisor * BigInteger(mid);
                //std::cout<<"mid: "<<mid<<std::endl;
                //std::cout<<"diff: "<<diff.to_string()<<std::endl;
                if (diff >= 0) {
                    count = mid;
                    //std::cout<<"count: "<<count<<std::endl;
                    left = mid + 1;
                } else {
                    right = mid - 1;
                }
                whileTest++;
                assert("operator/: while running too long! " && whileTest<10000);
            }
            quotient.digits[i] = count;
            current = current - divisor * BigInteger(count);
        }

        quotient.negative = (negative != other.negative);
        quotient.trim();
        
        //testing subdivision:
        //std::cout<<"dividend: "<<this->to_string()<<std::endl;
        //std::cout<<"divisor: "<<other.to_string()<<std::endl;
        //std::cout<<"quotient: "<<quotient.to_string()<<std::endl;
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
        if (digits.empty()) {
            oss << '0';
        } else {
            oss << digits.back();
            for (size_t i = digits.size() - 1; i-- > 0;) {
                oss << std::setw(DIGITS) << std::setfill('0') << digits[i];
            }
        }
        return oss.str();
    }

    friend std::ostream &operator<<(std::ostream &out, const BigInteger &value) {
        out << value.to_string();
        return out;
    }
};

BigInteger gcd(BigInteger a, BigInteger b) {
    int whileTest=0;
    while (b != BigInteger(0)) {
        BigInteger temp = b;
        b = a % b;
        a = temp;
        whileTest++;
        assert("gcd(): while running too long! " && whileTest<10000);
    }
    return a;
}
