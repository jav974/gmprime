#include <fstream>
#include <stdexcept>
#include <vector>
#include "gmprime.hpp"

using namespace GMPrime;

gmprime::gmprime(): _value(0) {
}

gmprime::gmprime(const std::string &value, const int base) {
    // Base 0 detects the base from the prefix (0x, 0b, 0), otherwise GMP supports bases 2 to 62
    if (base != 0 && (base < 2 || base > 62)) {
        throw std::invalid_argument("Invalid base " + std::to_string(base));
    }

    if (_value.set_str(value, base) != 0) {
        throw std::invalid_argument("Invalid number \"" + value + "\" in base " + std::to_string(base));
    }
}

mpz_class & gmprime::get_value() {
    return _value;
}

const mpz_class & gmprime::get_value() const {
    return _value;
}

namespace GMPrime {

gmprime operator+(const gmprime &lhs, const gmprime &rhs) {
    return gmprime(lhs) += rhs;
}

gmprime operator-(const gmprime &lhs, const gmprime &rhs) {
    return gmprime(lhs) -= rhs;
}

gmprime operator*(const gmprime &lhs, const gmprime &rhs) {
    return gmprime(lhs) *= rhs;
}

gmprime operator/(const gmprime &lhs, const gmprime &rhs) {
    return gmprime(lhs) /= rhs;
}

gmprime operator%(const gmprime &lhs, const gmprime &rhs) {
    return gmprime(lhs) %= rhs;
}

gmprime operator&(const gmprime &lhs, const gmprime &rhs) {
    return gmprime(lhs) &= rhs;
}

gmprime operator|(const gmprime &lhs, const gmprime &rhs) {
    return gmprime(lhs) |= rhs;
}

gmprime operator^(const gmprime &lhs, const gmprime &rhs) {
    return gmprime(lhs) ^= rhs;
}

bool operator==(const gmprime &lhs, const gmprime &rhs) {
    return lhs._value == rhs._value;
}

bool operator!=(const gmprime &lhs, const gmprime &rhs) {
    return lhs._value != rhs._value;
}

bool operator<=(const gmprime &lhs, const gmprime &rhs) {
    return lhs._value <= rhs._value;
}

bool operator>=(const gmprime &lhs, const gmprime &rhs) {
    return lhs._value >= rhs._value;
}

bool operator<(const gmprime &lhs, const gmprime &rhs) {
    return lhs._value < rhs._value;
}

bool operator>(const gmprime &lhs, const gmprime &rhs) {
    return lhs._value > rhs._value;
}

} // GMPrime

gmprime & gmprime::operator+=(const gmprime &other) {
    _value += other._value;
    return *this;
}

gmprime & gmprime::operator-=(const gmprime &other) {
    _value -= other._value;
    return *this;
}

gmprime & gmprime::operator*=(const gmprime &other) {
    _value *= other._value;
    return *this;
}

// GMP raises SIGFPE on division by zero
static void throw_if_zero_divisor(const mpz_class &divisor) {
    if (mpz_sgn(divisor.get_mpz_t()) == 0) {
        throw std::domain_error("Division by zero");
    }
}

gmprime & gmprime::operator/=(const gmprime &other) {
    throw_if_zero_divisor(other._value);
    _value /= other._value;
    return *this;
}

gmprime & gmprime::operator%=(const gmprime &other) {
    throw_if_zero_divisor(other._value);
    _value %= other._value;
    return *this;
}

gmprime & gmprime::operator++() {
    ++_value;
    return *this;
}

gmprime & gmprime::operator--() {
    --_value;
    return *this;
}

gmprime & gmprime::operator&=(const gmprime &other) {
    _value &= other._value;
    return *this;
}

gmprime & gmprime::operator|=(const gmprime &other) {
    _value |= other._value;
    return *this;
}

gmprime & gmprime::operator^=(const gmprime &other) {
    _value ^= other._value;
    return *this;
}

gmprime & gmprime::operator<<=(const size_t shift) {
    _value <<= shift;
    return *this;
}

gmprime & gmprime::operator>>=(const size_t shift) {
    _value >>= shift;
    return *this;
}

gmprime gmprime::operator<<(const size_t shift) const {
    return gmprime(*this) <<= shift;
}

gmprime gmprime::operator>>(const size_t shift) const {
    return gmprime(*this) >>= shift;
}

gmprime gmprime::operator~() const {
    gmprime result(*this);
    result._value = ~result._value;
    return result;
}

gmprime gmprime::operator+() const {
    gmprime result(*this);
    result._value = +_value;
    return result;
}

gmprime gmprime::operator-() const {
    gmprime result(*this);
    result._value = -result._value;
    return result;
}

size_t gmprime::used_bits() const {
    if (_value == 0) return 0;

    return size_in_base(2);
}

size_t gmprime::size_in_base(const int base) const {
    return mpz_sizeinbase(_value.get_mpz_t(), base);
}

bool gmprime::get_bit(const size_t pos) const {
    // mpz_tstbit uses two's complement for negative values, read the magnitude instead
    const auto limb = mpz_getlimbn(_value.get_mpz_t(), static_cast<mp_size_t>(pos / GMP_NUMB_BITS));
    return (limb >> (pos % GMP_NUMB_BITS)) & 1;
}

gmprime & gmprime::set_bit(const size_t pos, const bool value) {
    // mpz_setbit/mpz_clrbit use two's complement for negative values: work on the magnitude (mpz_neg is O(1))
    const bool negative = mpz_sgn(_value.get_mpz_t()) < 0;

    if (negative) {
        mpz_neg(_value.get_mpz_t(), _value.get_mpz_t());
    }

    if (value) {
        mpz_setbit(_value.get_mpz_t(), pos);
    } else {
        mpz_clrbit(_value.get_mpz_t(), pos);
    }

    if (negative) {
        mpz_neg(_value.get_mpz_t(), _value.get_mpz_t());
    }

    return *this;
}

std::string gmprime::to_string(const int base) const {
    return _value.get_str(base);
}

gmprime & gmprime::load_from_file(const std::string &file_path) {
    // Open the file in binary mode
    std::ifstream file(file_path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Error opening file");
    }

    // Determine the size of the file
    file.seekg(0, std::ios::end);
    const std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    // Read the file into a buffer
    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size)) {
        throw std::runtime_error("Error reading file");
    }

    file.close();

    // Import the buffer into the mpz_t variable
    mpz_import(_value.get_mpz_t(), size / sizeof(char), -1, sizeof(char), 0, 0, buffer.data());
    _byte_size = static_cast<size_t>(size);

    return *this;
}

void gmprime::save_to_file(const std::string &file_path) const {
    // Open the file in binary mode
    std::ofstream file(file_path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Error opening file");
    }

    // Determine the size of the export buffer
    size_t size;
    void* data = mpz_export(nullptr, &size, -1, 1, 1, 0, _value.get_mpz_t());

    // Write the data to the file
    file.write(static_cast<const char*>(data), static_cast<std::streamsize>(size));

    if (_byte_size > size) {
        const std::vector<char> padding(_byte_size - size, 0);
        file.write(padding.data(), static_cast<std::streamsize>(padding.size()));
    }
    file.close();

    // Free the allocated memory for export
    free(data);
}

gmprime gmprime::sqrt() const {
    // GMP raises SIGFPE on a negative value
    if (mpz_sgn(_value.get_mpz_t()) < 0) {
        throw std::domain_error("Square root of a negative number");
    }

    gmprime result;
    mpz_sqrt(result._value.get_mpz_t(), _value.get_mpz_t());
    return result;
}

gmprime gmprime::pow(const unsigned long int n) const {
    gmprime result;
    mpz_pow_ui(result._value.get_mpz_t(), _value.get_mpz_t(), n);
    return result;
}

bool gmprime::is_perfect_power() const {
    return mpz_perfect_power_p(_value.get_mpz_t()) != 0;
}

bool gmprime::is_divisible_by_pow2(const unsigned long int n) const {
    return mpz_divisible_2exp_p(_value.get_mpz_t(), n) != 0;
}

bool gmprime::is_divisible_by(const gmprime &other) const {
    return mpz_divisible_p(_value.get_mpz_t(), other._value.get_mpz_t()) != 0;
}

gmprime gmprime::next_prime() const {
    gmprime result(*this);
    mpz_nextprime(result._value.get_mpz_t(), _value.get_mpz_t());
    return result;
}

gmprime gmprime::factorial(const unsigned int n) {
    gmprime result;
    mpz_fac_ui(result._value.get_mpz_t(), n);
    return result;
}
