#include <fstream>
#include <vector>
#include "gmprime.hpp"

using namespace GMPrime;

gmprime::gmprime(const std::string &value, const int base) {
    _value.set_str(value, base);
}

mpz_class & gmprime::get_value() {
    return _value;
}

const mpz_class & gmprime::get_value() const {
    return _value;
}

gmprime gmprime::operator+(const gmprime &other) const {
    return gmprime(*this) += other;
}

gmprime gmprime::operator-(const gmprime &other) const {
    return gmprime(*this) -= other;
}

gmprime gmprime::operator*(const gmprime &other) const {
    return gmprime(*this) *= other;
}

gmprime gmprime::operator/(const gmprime &other) const {
    return gmprime(*this) /= other;
}

gmprime gmprime::operator%(const gmprime &other) const {
    return gmprime(*this) %= other;
}

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

gmprime & gmprime::operator/=(const gmprime &other) {
    _value /= other._value;
    return *this;
}

gmprime & gmprime::operator%=(const gmprime &other) {
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

gmprime gmprime::operator&(const gmprime &other) const {
    return gmprime(*this) &= other;
}

gmprime gmprime::operator|(const gmprime &other) const {
    return gmprime(*this) |= other;
}

gmprime gmprime::operator^(const gmprime &other) const {
    return gmprime(*this) ^= other;
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

bool gmprime::operator==(const gmprime &other) const {
    return _value == other._value;
}

bool gmprime::operator!=(const gmprime &other) const {
    return _value != other._value;
}

bool gmprime::operator<=(const gmprime &other) const {
    return _value <= other._value;
}

bool gmprime::operator>=(const gmprime &other) const {
    return _value >= other._value;
}

bool gmprime::operator<(const gmprime &other) const {
    return _value < other._value;
}

bool gmprime::operator>(const gmprime &other) const {
    return _value > other._value;
}

size_t gmprime::used_bits() const {
    if (_value == 0) return 0;

    return size_in_base(2);
}

size_t gmprime::size_in_base(const int base) const {
    return mpz_sizeinbase(_value.get_mpz_t(), base);
}

bool gmprime::get_bit(const size_t pos) const {
    return mpz_tstbit(_value.get_mpz_t(), pos) == 1;
}

gmprime & gmprime::set_bit(const size_t pos, const bool value) {
    if (value) {
        mpz_setbit(_value.get_mpz_t(), pos);
    } else {
        mpz_clrbit(_value.get_mpz_t(), pos);
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
    file.close();

    // Free the allocated memory for export
    free(data);
}

gmprime gmprime::sqrt() const {
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
