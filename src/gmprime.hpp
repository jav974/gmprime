#ifndef GMPRIME_HPP
#define GMPRIME_HPP

#include <functional>
#include <cstdint>
#include <gmpxx.h>
#include "gmprime_ios.hpp"
#include "gmprime_literals.hpp"

namespace GMPrime {

template<typename T = mp_limb_t, std::enable_if_t<std::is_integral_v<T>, bool> = true>
class reference;

template<typename T = mp_limb_t, std::enable_if_t<std::is_integral_v<T>, bool> = true>
class word_adapter;

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

    explicit gmprime(const std::string& value, int base = 10);

    mpz_class& get_value();
    [[nodiscard]] const mpz_class& get_value() const;

    // ASSIGNMENT OPERATORS =
    gmprime& operator=(const gmprime& other) = default;
    gmprime& operator=(gmprime&& other) noexcept = default;

    // ARITHMETIC OPERATORS +=  -=  *=  /=  %=  +   -   *   /   %   ++  --
    gmprime operator+(const gmprime& other) const;
    gmprime operator-(const gmprime& other) const;
    gmprime operator*(const gmprime& other) const;
    gmprime operator/(const gmprime& other) const;
    gmprime operator%(const gmprime& other) const;
    gmprime& operator+=(const gmprime& other);
    gmprime& operator-=(const gmprime& other);
    gmprime& operator*=(const gmprime& other);
    gmprime& operator/=(const gmprime& other);
    gmprime& operator%=(const gmprime& other);
    gmprime& operator++();
    gmprime& operator--();

    // BINARY OPERATORS &=  &   |=  |   ^=  ^   <<= <<  >>= >>  ~
    gmprime& operator&=(const gmprime& other);
    gmprime& operator|=(const gmprime& other);
    gmprime& operator^=(const gmprime& other);
    gmprime& operator<<=(size_t shift);
    gmprime& operator>>=(size_t shift);
    gmprime operator&(const gmprime& other) const;
    gmprime operator|(const gmprime& other) const;
    gmprime operator^(const gmprime& other) const;
    gmprime operator<<(size_t shift) const;
    gmprime operator>>(size_t shift) const;
    gmprime operator~() const;

    // UNARY OPERATORS -    +
    gmprime operator+() const;
    gmprime operator-() const;

    // COMPARISON OPERATORS  ==  <=  >=  <  >  !=
    bool operator==(const gmprime &other) const;
    bool operator!=(const gmprime &other) const;
    bool operator<=(const gmprime &other) const;
    bool operator>=(const gmprime &other) const;
    bool operator<(const gmprime &other) const;
    bool operator>(const gmprime &other) const;

    template<typename word_t = mp_limb_t, std::enable_if_t<std::is_integral_v<word_t>, bool> = true>
    word_t operator[](long int index) const;

    template<typename word_t>
    word_t operator[](const word_index<word_t>& index) const;

    template<typename word_t = mp_limb_t, std::enable_if_t<std::is_integral_v<word_t>, bool> = true>
    reference<word_t> operator[](long int index);

    template<typename word_t = mp_limb_t, std::enable_if_t<std::is_integral_v<word_t>, bool> = true>
    word_adapter<word_t> as_span(bool reverse = false);

    template<typename word_t = unsigned long int>
    [[nodiscard]] size_t used_words() const;

    [[nodiscard]] size_t used_bits() const;
    [[nodiscard]] size_t size_in_base(int base = 2) const;
    [[nodiscard]] bool get_bit(size_t pos) const;
    gmprime& set_bit(size_t pos, bool value);

    template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    void set_partial_limb(uint64_t index, word_t value);

    template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    word_t get_partial_limb(uint64_t index) const;

    [[nodiscard]] std::string to_string(int base = 10) const;

    template<typename word_t, std::enable_if_t<std::is_arithmetic_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    word_t to_primitive() const;

    gmprime& load_from_file(const std::string& file_path);
    void save_to_file(const std::string& file_path) const;

    [[nodiscard]] gmprime sqrt() const;
    [[nodiscard]] gmprime pow(unsigned long int n) const;

    [[nodiscard]] bool is_perfect_power() const;
    [[nodiscard]] bool is_divisible_by_pow2(unsigned long int n) const;
    [[nodiscard]] bool is_divisible_by(const gmprime& other) const;
    [[nodiscard]] gmprime next_prime() const;
    [[nodiscard]] static gmprime factorial(unsigned int n);
private:
    mpz_class _value;
};

#include "gmprime_reference.tpp"
#include "gmprime_adapter.tpp"
#include "gmprime.tpp"

} // GMPrime

#include "gmprime_stream.tpp" // NOLINT

#endif //GMPRIME_HPP
