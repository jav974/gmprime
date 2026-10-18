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

void gmprime::shift_left(const uintmax_t count) {
    const auto z = _value.get_mpz_t();

    if (mpz_sgn(z) == 0) {
        return;
    }

    // GMP aborts the program beyond its maximum number of limbs (same limit as in its _mpz_realloc), and
    // takes shift counts as unsigned long (32 bits on Windows)
    constexpr uintmax_t max_limbs = sizeof(mp_size_t) == sizeof(int)
        ? std::numeric_limits<unsigned long>::max() / GMP_NUMB_BITS
        : static_cast<uintmax_t>(std::numeric_limits<int>::max());
    constexpr uintmax_t max_bits = max_limbs * GMP_NUMB_BITS;

    if (count > std::numeric_limits<mp_bitcnt_t>::max() || count > max_bits - used_bits()) {
        throw std::length_error("Left shift beyond the maximum size of a number");
    }

    mpz_mul_2exp(z, z, static_cast<mp_bitcnt_t>(count));
}

void gmprime::shift_right(const uintmax_t count) {
    const auto z = _value.get_mpz_t();

    // Every bit shifted out: 0, or -1 for a negative number since the shift rounds toward minus infinity.
    // Also keeps a huge count from being truncated to unsigned long.
    if (count >= used_bits()) {
        mpz_set_si(z, mpz_sgn(z) < 0 ? -1 : 0);
        return;
    }

    mpz_fdiv_q_2exp(z, z, static_cast<mp_bitcnt_t>(count));
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

bit_window<false> gmprime::bits(const size_t pos, const size_t width) {
    return {*this, pos, width};
}

bit_window<true> gmprime::bits(const size_t pos, const size_t width) const {
    return {*this, pos, width};
}

gmprime gmprime::read_window(const size_t pos, const size_t width) const {
    gmprime result;
    const auto r = result._value.get_mpz_t();

    // Truncating division shifts the magnitude, whatever the sign, which is then dropped
    mpz_tdiv_q_2exp(r, _value.get_mpz_t(), pos);
    mpz_abs(r, r);
    mpz_fdiv_r_2exp(r, r, width);

    return result;
}

void gmprime::write_window(const size_t pos, const size_t width, const gmprime &value) {
    // Rebuilds the magnitude from its bits below and above the window and the low width bits of value
    // in between, then restores the sign. value is read first: it may be this number itself.
    const auto z = _value.get_mpz_t();
    const bool negative = mpz_sgn(z) < 0;
    mpz_class middle, low, high;

    mpz_abs(middle.get_mpz_t(), value._value.get_mpz_t());
    mpz_fdiv_r_2exp(middle.get_mpz_t(), middle.get_mpz_t(), width);

    mpz_abs(z, z);
    mpz_fdiv_r_2exp(low.get_mpz_t(), z, pos);
    mpz_fdiv_q_2exp(high.get_mpz_t(), z, pos + width);

    mpz_mul_2exp(high.get_mpz_t(), high.get_mpz_t(), width);
    mpz_ior(high.get_mpz_t(), high.get_mpz_t(), middle.get_mpz_t());
    mpz_mul_2exp(high.get_mpz_t(), high.get_mpz_t(), pos);
    mpz_ior(z, high.get_mpz_t(), low.get_mpz_t());

    if (negative) {
        mpz_neg(z, z);
    }
}

namespace {

// mpn_popcount and mpn_hamdist return an unsigned long, 32 bits on Windows: count by chunks small enough
constexpr auto max_counted_limbs = static_cast<size_t>(std::min<uintmax_t>(
    std::numeric_limits<mp_bitcnt_t>::max() / GMP_NUMB_BITS,
    static_cast<uintmax_t>(std::numeric_limits<mp_size_t>::max())));

size_t popcount_limbs(const mp_limb_t *limbs, size_t count) {
    size_t total = 0;

    while (count > 0) {
        const auto chunk = std::min(count, max_counted_limbs);
        total += mpn_popcount(limbs, static_cast<mp_size_t>(chunk));
        limbs += chunk;
        count -= chunk;
    }

    return total;
}

size_t hamdist_limbs(const mp_limb_t *a, const mp_limb_t *b, size_t count) {
    size_t total = 0;

    while (count > 0) {
        const auto chunk = std::min(count, max_counted_limbs);
        total += mpn_hamdist(a, b, static_cast<mp_size_t>(chunk));
        a += chunk;
        b += chunk;
        count -= chunk;
    }

    return total;
}

// Index of the lowest bit set in a non-zero limb
size_t lowest_set_bit(mp_limb_t limb) {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<size_t>(__builtin_ctzll(limb));
#else
    size_t index = 0;
    while ((limb & 1) == 0) {
        limb >>= 1;
        ++index;
    }
    return index;
#endif
}

} // namespace

// mpz_popcount, mpz_hamdist, mpz_scan0 and mpz_scan1 use two's complement for negative values (an infinity of
// bits set), so these work on the limbs of the magnitude directly

size_t gmprime::popcount() const {
    const auto z = _value.get_mpz_t();
    return popcount_limbs(mpz_limbs_read(z), mpz_size(z));
}

size_t gmprime::hamming_distance(const gmprime &other) const {
    const auto a = _value.get_mpz_t();
    const auto b = other._value.get_mpz_t();
    const auto size_a = mpz_size(a);
    const auto size_b = mpz_size(b);
    const auto common = std::min(size_a, size_b);
    const auto longer = size_a > size_b ? mpz_limbs_read(a) : mpz_limbs_read(b);

    // Beyond the shorter magnitude, every bit set in the longer one differs
    return hamdist_limbs(mpz_limbs_read(a), mpz_limbs_read(b), common)
           + popcount_limbs(longer + common, std::max(size_a, size_b) - common);
}

size_t gmprime::next_set_bit(const size_t from) const {
    constexpr size_t limb_bits = GMP_NUMB_BITS;
    const auto z = _value.get_mpz_t();
    const auto size = mpz_size(z);
    const auto limbs = mpz_limbs_read(z);
    auto index = from / limb_bits;

    if (index >= size) {
        return npos;
    }

    // Ignores the bits below from in its limb, then skips zero limbs whole
    auto limb = limbs[index] & (~static_cast<mp_limb_t>(0) << (from % limb_bits));

    while (limb == 0) {
        if (++index == size) {
            return npos;
        }
        limb = limbs[index];
    }

    return index * limb_bits + lowest_set_bit(limb);
}

size_t gmprime::next_clear_bit(const size_t from) const {
    constexpr size_t limb_bits = GMP_NUMB_BITS;
    const auto z = _value.get_mpz_t();
    const auto size = mpz_size(z);
    const auto limbs = mpz_limbs_read(z);
    auto index = from / limb_bits;

    if (index >= size) {
        return from;
    }

    // Same as next_set_bit on the inverted limbs, the bits beyond the magnitude being clear
    auto limb = ~limbs[index] & (~static_cast<mp_limb_t>(0) << (from % limb_bits));

    while (limb == 0) {
        if (++index == size) {
            return size * limb_bits;
        }
        limb = ~limbs[index];
    }

    return index * limb_bits + lowest_set_bit(limb);
}

set_bit_range gmprime::set_bits() const & {
    return set_bit_range(*this);
}

set_bit_range gmprime::set_bits() && {
    return set_bit_range(std::move(*this));
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

gmprime GMPrime::abs(const gmprime &value) {
    gmprime result;
    mpz_abs(result.get_value().get_mpz_t(), value.get_value().get_mpz_t());
    return result;
}

gmprime GMPrime::gcd(const gmprime &a, const gmprime &b) {
    gmprime result;
    mpz_gcd(result.get_value().get_mpz_t(), a.get_value().get_mpz_t(), b.get_value().get_mpz_t());
    return result;
}

gmprime GMPrime::lcm(const gmprime &a, const gmprime &b) {
    gmprime result;
    mpz_lcm(result.get_value().get_mpz_t(), a.get_value().get_mpz_t(), b.get_value().get_mpz_t());
    return result;
}

std::string GMPrime::to_string(const gmprime &value) {
    return value.to_string();
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
