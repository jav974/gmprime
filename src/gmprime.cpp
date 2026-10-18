#include <filesystem>
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

gmprime::operator bool() const {
    return mpz_sgn(_value.get_mpz_t()) != 0;
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

    // Determine the size of the file. Not with seekg/tellg: they give a meaningless size for anything
    // but a regular file (e.g. on Linux, opening a directory succeeds and tellg() returns LLONG_MAX on ext4)
    std::error_code error;
    const auto size = static_cast<std::streamsize>(std::filesystem::file_size(file_path, error));
    if (error) {
        throw std::runtime_error("Error reading file size: " + error.message());
    }

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

    // Export into our own buffer: a buffer allocated by GMP would have to be released through
    // GMP's memory functions, which may be customized (mp_set_memory_functions)
    const auto mpz = _value.get_mpz_t();
    std::vector<char> buffer(mpz_sgn(mpz) == 0 ? 0 : (mpz_sizeinbase(mpz, 2) + CHAR_BIT - 1) / CHAR_BIT);
    size_t size = 0;
    if (!buffer.empty()) {
        mpz_export(buffer.data(), &size, -1, 1, 1, 0, mpz);
    }

    // Write the data to the file
    file.write(buffer.data(), static_cast<std::streamsize>(size));

    if (_byte_size > size) {
        const std::vector<char> padding(_byte_size - size, 0);
        file.write(padding.data(), static_cast<std::streamsize>(padding.size()));
    }

    file.close();
    if (!file) {
        throw std::runtime_error("Error writing file");
    }
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

size_t std::hash<gmprime>::operator()(const gmprime &value) const noexcept {
    // GMP keeps numbers normalized (no high zero limb), so equal values have the same limbs
    const auto mpz = value.get_value().get_mpz_t();
    const auto limbs = mpz_limbs_read(mpz);
    const auto size = mpz_size(mpz);
    auto seed = std::hash<int>{}(mpz_sgn(mpz));

    // Same mixing as boost::hash_combine
    for (size_t i = 0; i < size; ++i) {
        seed ^= std::hash<mp_limb_t>{}(limbs[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    return seed;
}
