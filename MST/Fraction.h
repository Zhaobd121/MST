#pragma once
#include <cstddef>    
#include <stdexcept>  
#include <ostream>

class Fraction {
public:

    constexpr Fraction(int n = 0, int d = 1) : num_(n), den_(d) {
        if (den_ == 0) throw std::invalid_argument("denominator must be > 0");
    }

    constexpr int num() const noexcept { return num_; }
    constexpr int den() const noexcept { return den_; }

    friend constexpr bool operator==(const Fraction& a, const Fraction& b) noexcept {
        return static_cast<long long>(a.num_) * b.den_
            == static_cast<long long>(b.num_) * a.den_;
    }
    friend constexpr bool operator!=(const Fraction& a, const Fraction& b) noexcept {
        return !(a == b);
    }

    friend constexpr bool operator<(const Fraction& a, const Fraction& b) noexcept {
        return static_cast<long long>(a.num_) * b.den_
             < static_cast<long long>(b.num_) * a.den_;
    }
    friend constexpr bool operator>(const Fraction& a, const Fraction& b) noexcept { return b < a; }
    friend constexpr bool operator<=(const Fraction& a, const Fraction& b) noexcept { return !(b < a); }
    friend constexpr bool operator>=(const Fraction& a, const Fraction& b) noexcept { return !(a < b); }

    friend std::ostream& operator<<(std::ostream& os, const Fraction& f) {
        if (f.den_ == 1)
            os << f.num_;                
        else
            os << f.num_ << '/' << f.den_; 
        return os;
    }

private:
    int num_, den_;
};