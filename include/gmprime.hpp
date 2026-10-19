#ifndef GMPRIME_HPP
#define GMPRIME_HPP

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
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

class set_bit_range;

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
    // Shift counts of any integer type. A negative count throws std::domain_error, a left shift beyond the
    // maximum size of a GMP number throws std::length_error (GMP would abort the program). Right shifts round
    // toward minus infinity, like an int.
    template<typename shift_t, std::enable_if_t<std::is_integral_v<shift_t> && !std::is_same_v<shift_t, bool>, bool> = true>
    gmprime& operator<<=(shift_t shift);
    template<typename shift_t, std::enable_if_t<std::is_integral_v<shift_t> && !std::is_same_v<shift_t, bool>, bool> = true>
    gmprime& operator>>=(shift_t shift);
    friend gmprime operator&(const gmprime& lhs, const gmprime& rhs);
    friend gmprime operator|(const gmprime& lhs, const gmprime& rhs);
    friend gmprime operator^(const gmprime& lhs, const gmprime& rhs);
    template<typename shift_t, std::enable_if_t<std::is_integral_v<shift_t> && !std::is_same_v<shift_t, bool>, bool> = true>
    gmprime operator<<(shift_t shift) const;
    template<typename shift_t, std::enable_if_t<std::is_integral_v<shift_t> && !std::is_same_v<shift_t, bool>, bool> = true>
    gmprime operator>>(shift_t shift) const;
    gmprime operator~() const;

    // UNARY OPERATORS -    +
    gmprime operator+() const;
    gmprime operator-() const;

    // CONVERSION OPERATOR bool: true when non-zero. Explicit so that if (a) and !a work like with
    // a primitive, without a gmprime silently converting to bool in arithmetic or comparisons
    explicit operator bool() const;

    // COMPARISON OPERATORS  ==  <=  >=  <  >  !=
    // With a floating point value, the comparison is exact, without converting either side (see gmprime.tpp):
    // gmprime(5) < 5.5 is true, and NaN compares unordered (every comparison false, except !=)
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
    // Bit access works on the magnitude, as GMP stores it, consistently with word access ([], as_span, bits).
    // Bitwise operators (& | ^ ~ >>) follow two's complement instead, like an int: for a negative number,
    // get_bit(k) can differ from ((a >> k) & 1), e.g. for -2, bit 2 is 0 in the magnitude and 1 in two's complement.
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

    // BIT ANALYSIS, on the magnitude like every bit access
    static constexpr size_t npos = std::numeric_limits<size_t>::max();
    // Number of bits set
    [[nodiscard]] size_t popcount() const;
    // Number of bits that differ from other (Hamming distance)
    [[nodiscard]] size_t hamming_distance(const gmprime& other) const;
    // Position of the first bit set at or after from, npos if there is none
    [[nodiscard]] size_t next_set_bit(size_t from = 0) const;
    // Position of the first clear bit at or after from: always found, the magnitude being followed by zeros
    [[nodiscard]] size_t next_clear_bit(size_t from = 0) const;
    // Positions of the bits set, in increasing order: for (const size_t pos : a.set_bits()) ...
    // Skips clear bits a limb at a time, so its cost follows the number of bits set more than the size.
    // Iterators are invalidated by any change to the number. On a temporary, the range keeps the number,
    // so that for (const size_t pos : (a & mask).set_bits()) does not dangle.
    [[nodiscard]] set_bit_range set_bits() const &;
    [[nodiscard]] set_bit_range set_bits() &&;

    // Occurrences of each symbol of width bits (1 to max_histogram_width): result[s] counts the windows
    // holding the value s, bit i of a symbol starting at bit start being bit start + i of the number (so a
    // file's bytes read as little-endian). Windows start every step bits from bit 0: step = width (default)
    // for aligned symbols, e.g. bytes with width = 8, or step = 1 for a sliding window. Only complete windows
    // within length bits count. length defaults to the whole number in whole bytes, like a file (the high byte
    // of a text is partial, ASCII leaving its top bit clear), and to the size of the last loaded file if larger:
    // its high zero bytes are part of it.
    // Throws std::invalid_argument for an invalid width.
    static constexpr size_t max_histogram_width = 24;
    [[nodiscard]] std::vector<size_t> histogram(size_t width, size_t step = 0, size_t length = 0) const;
    // Shannon entropy of a histogram, in bits per symbol: a lower bound on the average size of a symbol
    // for any code that encodes symbols independently. 0 for an empty histogram.
    [[nodiscard]] static double entropy(const std::vector<size_t>& histogram);

    // CONTEXT ANALYSIS, on aligned symbols of width bits (1 to 32), over the same bits as histogram():
    // bytes of a text with width = 8. Throws std::invalid_argument for an invalid width.
    struct context_entropy {
        double naive;      // bits per symbol, from the observed frequencies: optimistic when contexts are rare
        double corrected;  // bits per symbol, with the Miller-Madow correction of part of that bias
        size_t contexts;   // distinct contexts seen
        size_t samples;    // symbols measured: those with order symbols before them
    };
    // Entropy of a symbol given the order symbols before it; order 0 gives the entropy of histogram(width)
    [[nodiscard]] context_entropy conditional_entropy(size_t order, size_t width = 8) const;
    // Bits an adaptive order-k context model takes to code every symbol, learning the frequencies of each
    // context as it goes (Krichevsky-Trofimov estimator; the first order symbols cost width bits each).
    // Unlike an entropy estimate, an achievable size, learning included: comparable to real compressors.
    [[nodiscard]] double adaptive_code_length(size_t order, size_t width = 8) const;

    // REPETITIONS, on the bytes of the number from the lowest (a loaded file's bytes, in order)
    // Longest previous factor: result[i] is the length of the longest string starting at byte i that also
    // starts at an earlier byte (they may overlap), i.e. what an LZ77 compressor with an unbounded window
    // can reference. Linear time from a suffix array; about 20 bytes of memory per byte.
    [[nodiscard]] std::vector<uint32_t> longest_previous_factors() const;
    struct repetition_stats {
        size_t bytes;          // bytes analyzed
        size_t matches;        // repeats referenced
        size_t matched_bytes;  // bytes covered by them
        size_t literals;       // bytes left: bytes - matched_bytes
    };
    // Greedy LZ77 parse: at each byte, a match if a repeat of at least min_length bytes starts there
    // (the longest one), else a literal
    [[nodiscard]] repetition_stats repetitions(size_t min_length = 4) const;

    // FOLDS: self-similarity of the magnitude at every distance up to max_distance bits. result[d] counts
    // the bits p (over the analyzed bits, p + d included) equal to bit p + d; result[0] is the number of bits.
    // Peaks of result[d] / (result[0] - d) reveal fold widths: a bitmap's row, a record, a line. Computed on
    // the whole number with one shift, XOR and popcount per distance (bit order: gmprime's, each byte from
    // its lowest bit).
    [[nodiscard]] std::vector<size_t> autocorrelation(size_t max_distance) const;

    // The same bytes with the bits of each byte reversed: bit p of the result is bit p of the bytes read as a
    // stream, each byte from its highest bit, the order of a bitmap's pixels (and of coding::stream_bit).
    // Shifts and autocorrelation of the result then work in stream order, where folds of any width line up.
    [[nodiscard]] gmprime stream_order() const;

    struct fold_analysis {
        size_t width = 0;        // best fold width in bits, 0 if the data is too short to search
        double score = 1;        // bits differing at that width, relative to independent bits: lower is more alike
        double median = 1;       // median score over the widths searched
        size_t layer = 0;        // rows per layer of the best 3D view, 3 or more (0: none)
        double layer_score = 1;  // score at width * layer

        // Whether the fold stands out from the other widths: its score below ratio times the median. 0.4
        // separates bitmaps and fixed-size records from text, code and executables in the Canterbury corpus.
        [[nodiscard]] bool strong(const double ratio = 0.4) const {
            return width > 0 && score < ratio * median;
        }
    };
    // Searches the fold widths from min_width to max_width bits (at most a quarter of the data), in stream
    // order: the width where bits most resemble the bits one row before, then the best layer of rows. With
    // sample_bytes > 0, only the first sample_bytes bytes are analyzed: the cost grows with the size, and the
    // first 128 KB find the same folds as whole files in the Canterbury corpus.
    [[nodiscard]] fold_analysis analyze_folds(size_t min_width = 64, size_t max_width = 8192, size_t sample_bytes = 0) const;

    template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    void set_partial_limb(uint64_t index, word_t value);

    template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    word_t get_partial_limb(uint64_t index) const;

    [[nodiscard]] std::string to_string(int base = 10) const;

    template<typename word_t, std::enable_if_t<std::is_arithmetic_v<word_t> && !std::is_same_v<word_t, bool>, bool> = true>
    word_t to_primitive() const;

    // Bytes of the number from the lowest, in whole bytes and up to the size of the last loaded file:
    // the contents of the file save_to_file would write
    [[nodiscard]] std::vector<unsigned char> bytes() const;
    // The number holding these bytes, the first one lowest, keeping high zero bytes like a loaded file
    [[nodiscard]] static gmprime from_bytes(const std::vector<unsigned char>& bytes);

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
    // C(n, k), the number of ways to choose k items among n (0 if k > n)
    [[nodiscard]] static gmprime binomial(unsigned long n, unsigned long k);

    // ENUMERATIVE CODING. A bit field of length bits with ones = popcount() bits set is one of C(length, ones)
    // arrangements. combination_rank is its index among them (combinatorial number system: the i-th bit set
    // from the lowest, i from 1, at position p, adds C(p, i)), below C(length, ones): about log2 C(length, ones)
    // bits, the information left in the field once length and ones are known. from_combination_rank rebuilds
    // the field. Both take O(length) steps on numbers of that size, so quadratic time: rank long fields in
    // blocks (rank_blocks in gmprime_coding.hpp).
    // Throws std::invalid_argument if a bit is set at or beyond length
    [[nodiscard]] gmprime combination_rank(size_t length) const;
    // Throws std::invalid_argument if ones exceeds length or rank is negative or not below C(length, ones)
    [[nodiscard]] static gmprime from_combination_rank(const gmprime& rank, size_t length, size_t ones);
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

    // Shift count as an unsigned integer, std::domain_error if negative
    template<typename shift_t>
    static uintmax_t shift_count(shift_t shift);
    void shift_left(uintmax_t count);
    void shift_right(uintmax_t count);

    // Number of bits the analysis functions cover: the used bits rounded up to whole bytes, or a loaded
    // file's size if larger
    [[nodiscard]] size_t analyzed_bits() const;
    // Aligned symbols of width bits (1 to 32) over analyzed_bits(), complete ones only
    [[nodiscard]] std::vector<uint32_t> symbols(size_t width) const;

    template<bool is_const>
    friend class bit_window;

    mpz_class _value;
    size_t _byte_size = 0;
};

#include "../src/gmprime_reference.tpp"
#include "../src/gmprime_adapter.tpp"
#include "../src/gmprime_bit_window.tpp"
#include "../src/gmprime_set_bits.hpp"
#include "../src/gmprime.tpp"

// Free functions, found by argument-dependent lookup like the ones for an int: abs(a), gcd(a, b).
// std::to_string cannot be overloaded (adding overloads to std is undefined behavior), so call
// to_string unqualified, after using std::to_string if needed.
[[nodiscard]] gmprime abs(const gmprime& value);
// Always >= 0, and gcd(0, 0) = 0, like std::gcd
[[nodiscard]] gmprime gcd(const gmprime& a, const gmprime& b);
// Always >= 0, and 0 if either is 0, like std::lcm
[[nodiscard]] gmprime lcm(const gmprime& a, const gmprime& b);
[[nodiscard]] std::string to_string(const gmprime& value);

} // GMPrime

namespace std {

// Hashes the value only (sign and limbs), consistently with operator==: the size of the last
// loaded file is not part of it
template<>
struct hash<GMPrime::gmprime> {
    size_t operator()(const GMPrime::gmprime& value) const noexcept;
};

// An unbounded integer: the members that only make sense for bounded or floating point types are zero,
// false, or gmprime(), like the primary template
template<>
class numeric_limits<GMPrime::gmprime> {
public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = true;
    static constexpr bool is_exact = true;
    static constexpr bool has_infinity = false;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
#if __cplusplus < 202302L
    static constexpr float_denorm_style has_denorm = denorm_absent;
    static constexpr bool has_denorm_loss = false;
#endif
    static constexpr float_round_style round_style = round_toward_zero;
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = false;
    static constexpr bool is_modulo = false;
    static constexpr int digits = 0;
    static constexpr int digits10 = 0;
    static constexpr int max_digits10 = 0;
    static constexpr int radix = 2;
    static constexpr int min_exponent = 0;
    static constexpr int min_exponent10 = 0;
    static constexpr int max_exponent = 0;
    static constexpr int max_exponent10 = 0;
    // Division by zero throws an exception rather than trapping
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;

    static GMPrime::gmprime min() { return {}; }
    static GMPrime::gmprime max() { return {}; }
    static GMPrime::gmprime lowest() { return {}; }
    static GMPrime::gmprime epsilon() { return {}; }
    static GMPrime::gmprime round_error() { return {}; }
    static GMPrime::gmprime infinity() { return {}; }
    static GMPrime::gmprime quiet_NaN() { return {}; }
    static GMPrime::gmprime signaling_NaN() { return {}; }
    static GMPrime::gmprime denorm_min() { return {}; }
};

} // std

#include "../src/gmprime_stream.tpp" // NOLINT

#endif //GMPRIME_HPP
