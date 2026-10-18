#ifndef GMPRIME_HPP
#define GMPRIME_HPP

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <gmpxx.h>
#include "../src/gmprime_ios.hpp"
#include "../src/gmprime_literals.hpp"

namespace GMPrime {

template<typename T = mp_limb_t, std::enable_if_t<std::is_integral_v<T>, bool> = true>
class reference;

template<typename T = mp_limb_t, std::enable_if_t<std::is_integral_v<T>, bool> = true>
class word_adapter;

template<bool is_const>
class bit_window;

class gmprime {
public:
    static constexpr auto bin = ios_flags{ios_base::bin};
    static constexpr auto oct = ios_flags{ios_base::oct};
    static constexpr auto dec = ios_flags{ios_base::dec};
    static constexpr auto hex = ios_flags{ios_base::hex};

    gmprime();
    gmprime(const gmprime& other) = default;
    gmprime(gmprime&& other) noexcept = default;

    template<typename word_t = mp_limb_t, std::enable_if_t<std::is_arithmetic_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    gmprime(word_t value); // NOLINT(*-explicit-constructor)

    // Throws std::invalid_argument if value is not a valid number in the given base
    explicit gmprime(const std::string& value, int base = 10);

    mpz_class& get_value();
    [[nodiscard]] const mpz_class& get_value() const;

    // ASSIGNMENT OPERATORS =
    gmprime& operator=(const gmprime& other) = default;
    gmprime& operator=(gmprime&& other) noexcept = default;

    // ARITHMETIC OPERATORS +=  -=  *=  /=  %=  +   -   *   /   %   ++  --
    // Binary operators are non-member friends so that a primitive left-hand operand converts too (e.g. 5 + a)
    // Division and modulo by zero throw std::domain_error
    friend gmprime operator+(const gmprime& lhs, const gmprime& rhs);
    friend gmprime operator-(const gmprime& lhs, const gmprime& rhs);
    friend gmprime operator*(const gmprime& lhs, const gmprime& rhs);
    friend gmprime operator/(const gmprime& lhs, const gmprime& rhs);
    friend gmprime operator%(const gmprime& lhs, const gmprime& rhs);
    gmprime& operator+=(const gmprime& other);
    gmprime& operator-=(const gmprime& other);
    gmprime& operator*=(const gmprime& other);
    gmprime& operator/=(const gmprime& other);
    gmprime& operator%=(const gmprime& other);
    gmprime& operator++();
    gmprime& operator--();
    // Postfix forms are deleted on purpose: they would copy the whole number to return the old value
    gmprime operator++(int) = delete;
    gmprime operator--(int) = delete;

    // BINARY OPERATORS &=  &   |=  |   ^=  ^   <<= <<  >>= >>  ~
    gmprime& operator&=(const gmprime& other);
    gmprime& operator|=(const gmprime& other);
    gmprime& operator^=(const gmprime& other);
    gmprime& operator<<=(size_t shift);
    gmprime& operator>>=(size_t shift);
    friend gmprime operator&(const gmprime& lhs, const gmprime& rhs);
    friend gmprime operator|(const gmprime& lhs, const gmprime& rhs);
    friend gmprime operator^(const gmprime& lhs, const gmprime& rhs);
    gmprime operator<<(size_t shift) const;
    gmprime operator>>(size_t shift) const;
    gmprime operator~() const;

    // UNARY OPERATORS -    +
    gmprime operator+() const;
    gmprime operator-() const;

    // CONVERSION OPERATOR bool: true when non-zero. Explicit so that if (a) and !a work like with
    // a primitive, without a gmprime silently converting to bool in arithmetic or comparisons
    explicit operator bool() const;

    // COMPARISON OPERATORS  ==  <=  >=  <  >  !=
    friend bool operator==(const gmprime& lhs, const gmprime& rhs);
    friend bool operator!=(const gmprime& lhs, const gmprime& rhs);
    friend bool operator<=(const gmprime& lhs, const gmprime& rhs);
    friend bool operator>=(const gmprime& lhs, const gmprime& rhs);
    friend bool operator<(const gmprime& lhs, const gmprime& rhs);
    friend bool operator>(const gmprime& lhs, const gmprime& rhs);

    template<typename word_t = mp_limb_t, std::enable_if_t<std::is_integral_v<word_t>, bool> = true>
    word_t operator[](long int index) const;

    template<typename word_t>
    word_t operator[](const word_index<word_t>& index) const;

    template<typename word_t = mp_limb_t, std::enable_if_t<std::is_integral_v<word_t>, bool> = true>
    reference<word_t> operator[](long int index);

    template<typename word_t>
    reference<word_t> operator[](const word_index<word_t>& index);

    template<typename word_t = mp_limb_t, std::enable_if_t<std::is_integral_v<word_t>, bool> = true>
    word_adapter<word_t> as_span(bool reverse = false);

    template<typename word_t = unsigned long int>
    [[nodiscard]] size_t used_words() const;

    [[nodiscard]] size_t used_bits() const;
    [[nodiscard]] size_t size_in_base(int base = 2) const;
    // Bit access works on the magnitude, consistently with word access ([] and as_span)
    [[nodiscard]] bool get_bit(size_t pos) const;
    gmprime& set_bit(size_t pos, bool value);

    // Window of width bits starting at bit pos, on the magnitude: a sign has no meaning within a window,
    // so only unsigned integers and gmprime convert to and from it (signed integers, floating point,
    // bool and characters are refused at compile time).
    //   uint16_t x = a.bits(1000, 13);  // std::overflow_error if width exceeds the type, whatever the value
    //   gmprime  g = a.bits(1000, 13);  // any width, always >= 0
    //   a.bits(1000, 13) = 0x1ABCu;     // keeps the low width bits of the value, like a C bit field
    //   a.bits(1000, 13) = g;           // same with the magnitude of a gmprime, any width
    // Writing keeps the sign of the number. Aligned word access is the special case bits(i * w, w).
    bit_window<false> bits(size_t pos, size_t width);
    [[nodiscard]] bit_window<true> bits(size_t pos, size_t width) const;

    template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    void set_partial_limb(uint64_t index, word_t value);

    template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    word_t get_partial_limb(uint64_t index) const;

    [[nodiscard]] std::string to_string(int base = 10) const;

    template<typename word_t, std::enable_if_t<std::is_arithmetic_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    word_t to_primitive() const;

    gmprime& load_from_file(const std::string& file_path);
    // Writes the magnitude: a file is a raw payload of bits, the sign is not stored.
    // Pads with zero bytes up to the size of the last loaded file: leading zeros of a number
    // (i.e. trailing zero bytes of the file) are not stored by GMP and would otherwise be lost.
    // This size is kept through in-place changes and copies, and reset by assigning a new value.
    void save_to_file(const std::string& file_path) const;

    // Throws std::domain_error for a negative value
    [[nodiscard]] gmprime sqrt() const;
    [[nodiscard]] gmprime pow(unsigned long int n) const;

    [[nodiscard]] bool is_perfect_power() const;
    [[nodiscard]] bool is_divisible_by_pow2(unsigned long int n) const;
    [[nodiscard]] bool is_divisible_by(const gmprime& other) const;
    [[nodiscard]] gmprime next_prime() const;
    [[nodiscard]] static gmprime factorial(unsigned int n);
private:
    // Bit field access on the magnitude, limb by limb with masks and shifts: independent of the
    // machine endianness and of the word width relative to the limb width. width must not exceed the bits of unsigned_t
    template<typename unsigned_t>
    [[nodiscard]] unsigned_t read_bits(size_t pos, size_t width) const;

    template<typename unsigned_t>
    void write_bits(size_t pos, size_t width, unsigned_t value);

    // Same on windows of any width, through GMP operations
    [[nodiscard]] gmprime read_window(size_t pos, size_t width) const;
    void write_window(size_t pos, size_t width, const gmprime& value);

    template<bool is_const>
    friend class bit_window;

    mpz_class _value;
    size_t _byte_size = 0;
};

#include "../src/gmprime_reference.tpp"
#include "../src/gmprime_adapter.tpp"
#include "../src/gmprime_bit_window.tpp"
#include "../src/gmprime.tpp"

} // GMPrime

namespace std {

// Hashes the value only (sign and limbs), consistently with operator==: the size of the last
// loaded file is not part of it
template<>
struct hash<GMPrime::gmprime> {
    size_t operator()(const GMPrime::gmprime& value) const noexcept;
};

} // std

#include "../src/gmprime_stream.tpp" // NOLINT

#endif //GMPRIME_HPP
