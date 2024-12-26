#include <gmp.h>
#include <iostream>

class Rational {
private:
    mpq_t value; // GMP rational number type

public:
    // Constructor
    Rational() {
        mpq_init(value);
    }

    // Constructor to initialize from numerator and denominator
    Rational(unsigned long numerator, unsigned long denominator) {
        mpq_init(value);
        mpq_set_ui(value, numerator, denominator);
        mpq_canonicalize(value);
    }

    // Copy constructor
    Rational(const Rational& other) {
        mpq_init(value);
        mpq_set(value, other.value);
    }

    // Destructor
    ~Rational() {
        mpq_clear(value);
    }

    // Overloaded addition operator
    Rational operator+(const Rational& other) const {
        Rational result;
        mpq_add(result.value, value, other.value);
        return result;
    }

    // Overloaded subtraction operator
    Rational operator-(const Rational& other) const {
        Rational result;
        mpq_sub(result.value, value, other.value);
        return result;
    }

    // Overloaded multiplication operator
    Rational operator*(const Rational& other) const {
        Rational result;
        mpq_mul(result.value, value, other.value);
        return result;
    }

    // Overloaded division operator
    Rational operator/(const Rational& other) const {
        Rational result;
        mpq_div(result.value, value, other.value);
        return result;
    }

    // Overloaded greater than operator
    bool operator>(const Rational& other) const {
        return mpq_cmp(value, other.value) > 0;
    }

    // Overloaded equality operator
    bool operator==(const Rational& other) const {
        return mpq_cmp(value, other.value) == 0;
    }

    // Output stream operator for printing
    friend std::ostream& operator<<(std::ostream& os, const Rational& r) {
        char* str = mpq_get_str(nullptr, 10, r.value);
        os << str;
        free(str); // Free the string allocated by GMP
        return os;
    }
};