#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <ctime>
#include <boost/ut.hpp>
#include <gmprime.hpp>
#include <gmprime_coding.hpp>
#include <gmprime_records.hpp>

using namespace boost::ut;
using namespace GMPrime;

template<typename L, typename R, typename = void>
struct has_non_member_operators : std::false_type {};

template<typename L, typename R>
struct has_non_member_operators<L, R, std::void_t<
    decltype(std::declval<L>() + std::declval<R>()),
    decltype(std::declval<L>() - std::declval<R>()),
    decltype(std::declval<L>() * std::declval<R>()),
    decltype(std::declval<L>() / std::declval<R>()),
    decltype(std::declval<L>() % std::declval<R>()),
    decltype(std::declval<L>() & std::declval<R>()),
    decltype(std::declval<L>() | std::declval<R>()),
    decltype(std::declval<L>() ^ std::declval<R>()),
    decltype(std::declval<L>() == std::declval<R>()),
    decltype(std::declval<L>() < std::declval<R>())>> : std::true_type {};

template<typename T, typename = void>
struct has_postfix_increment : std::false_type {};

template<typename T>
struct has_postfix_increment<T, std::void_t<decltype(std::declval<T&>()++)>> : std::true_type {};

template<typename T, typename = void>
struct has_postfix_decrement : std::false_type {};

template<typename T>
struct has_postfix_decrement<T, std::void_t<decltype(std::declval<T&>()--)>> : std::true_type {};

// GMP memory functions that count live allocations, delegating to the default ones
void* (*gmp_default_alloc)(size_t);
void* (*gmp_default_realloc)(void*, size_t, size_t);
void (*gmp_default_free)(void*, size_t);
long gmp_live_allocations = 0;

void* counting_alloc(const size_t size) {
    ++gmp_live_allocations;
    return gmp_default_alloc(size);
}

void* counting_realloc(void* ptr, const size_t old_size, const size_t new_size) {
    return gmp_default_realloc(ptr, old_size, new_size);
}

void counting_free(void* ptr, const size_t size) {
    --gmp_live_allocations;
    gmp_default_free(ptr, size);
}

// Positions of the bits set, read one by one: reference for the bit analysis functions
void naive_set_bits(const gmprime& value, std::vector<size_t>& positions) {
    for (size_t pos = 0; pos < value.used_bits(); ++pos) {
        if (value.get_bit(pos)) {
            positions.push_back(pos);
        }
    }
}

// The bytes of text, the first one lowest, like a loaded file without high zero bytes
void from_bytes(const std::string& text, gmprime& value) {
    for (size_t i = 0; i < text.size(); ++i) {
        value.bits(8 * i, 8) = static_cast<unsigned char>(text[i]);
    }
}

// Longest previous factor of each byte, by comparing with every earlier position
void naive_longest_previous_factors(const std::string& text, std::vector<uint32_t>& lpf) {
    for (size_t i = 0; i < text.size(); ++i) {
        uint32_t best = 0;
        for (size_t j = 0; j < i; ++j) {
            uint32_t length = 0;
            while (i + length < text.size() && text[i + length] == text[j + length]) {
                ++length;
            }
            best = std::max(best, length);
        }
        lpf.push_back(best);
    }
}

// Deterministic pseudo-random text over the first alphabet letters
void random_text(const size_t length, const unsigned alphabet, uint32_t seed, std::string& text) {
    for (size_t i = 0; i < length; ++i) {
        seed = seed * 1664525u + 1013904223u;
        text += static_cast<char>('a' + (seed >> 16) % alphabet);
    }
}

template<typename L>
void check_primitive_lhs(const L lhs) {
    if constexpr (has_non_member_operators<L, gmprime>::value) {
        const gmprime a(3);
        expect(lhs + a == 10);
        expect(lhs - a == 4);
        expect(lhs * a == 21);
        expect(lhs / a == 2);
        expect(lhs % a == 1);
        expect((lhs & a) == 3);
        expect((lhs | a) == 7);
        expect((lhs ^ a) == 4);
        expect(lhs == gmprime(7));
        expect(lhs > a);
    } else {
        expect(false) << "gmprime operators do not accept a primitive left-hand operand";
    }
}

suite<"gmprime::ctor"> ctor = [] {
    "ctor_signed_integer"_test = [] {
        const gmprime v(-42);
        expect(v == -42);
    };

    "ctor_unsigned_integer"_test = [] {
        const gmprime v(static_cast<unsigned long int>(42));
        expect(v == 42);
    };

    "ctor_floating_point"_test = [] {
        const gmprime v(-42.5);
        expect(v == -42);
    };

    "ctor_string_default"_test = [] {
        const gmprime v("42");
        expect(v == 42);
    };

    "ctor_string_base_2"_test = [] {
        const gmprime v("11", 2);
        expect(v == 3);
    };

    "ctor_string_base_8"_test = [] {
        const gmprime v("42", 8);
        expect(v == 34);
    };

    "ctor_string_base_16"_test = [] {
        const gmprime v("2a", 16);
        expect(v == 42);
    };

    "ctor_string_invalid"_test = [] {
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const gmprime v("12x4"); }));
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const gmprime v(""); }));
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const gmprime v("2", 2); }));
    };

    "ctor_string_invalid_base"_test = [] {
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const gmprime v("1", 1); }));
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const gmprime v("1", 63); }));
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const gmprime v("1", -2); }));
        expect(gmprime("0x2a", 0) == 42) << "base 0 detects the base from the prefix";
    };

    "ctor_floating_point_non_finite"_test = [] {
        expect(throws<std::domain_error>([] { [[maybe_unused]] const gmprime v(std::numeric_limits<double>::infinity()); }));
        expect(throws<std::domain_error>([] { [[maybe_unused]] const gmprime v(std::numeric_limits<double>::quiet_NaN()); }));
    };

    "ctor_64_bit_extremes"_test = [] {
        expect(gmprime(std::numeric_limits<int64_t>::min()) == gmprime("-9223372036854775808"));
        expect(gmprime(std::numeric_limits<int64_t>::max()) == gmprime("9223372036854775807"));
        expect(gmprime(std::numeric_limits<uint64_t>::max()) == gmprime("18446744073709551615"));
    };

#ifdef __SIZEOF_INT128__
    // Wider than long on every platform: exercises the path used for 64-bit integers where long is 32 bits (e.g. Windows)
    "ctor_wider_than_long"_test = [] {
        const auto max = ~static_cast<unsigned __int128>(0);
        expect(gmprime(max) == gmprime("340282366920938463463374607431768211455"));

        const auto min = -static_cast<__int128>(max >> 1) - 1;
        expect(gmprime(min) == gmprime("-170141183460469231731687303715884105728"));
        expect(gmprime(static_cast<__int128>(-42)) == -42);
    };
#endif
};

suite<"gmprime::assigment"> assign = [] {
    "=arithmetic"_test = [] {
        const gmprime v = 23156;
        expect(v == 23156);
    };
};

suite<"gmprime::comparison"> comparison = [] {
    "<"_test = [] {
        const gmprime lhs(100);
        gmprime rhs(50);
        expect((lhs < rhs) == false);
        expect(rhs < lhs);
        rhs = 100;
        expect((lhs < rhs) == false);
    };

    "<="_test = [] {
        gmprime lhs(100), rhs(50);
        expect((lhs <= rhs) == false);
        expect(rhs <= lhs);
        rhs = 100;
        expect(lhs <= rhs);
    };

    ">"_test = [] {
        gmprime lhs(100), rhs(50);
        expect(lhs > rhs);
        expect((rhs > lhs) == false);
        rhs = 100;
        expect((lhs > rhs) == false);
    };

    ">="_test = [] {
        gmprime lhs(100), rhs(50);
        expect(lhs >= rhs);
        expect((rhs >= lhs) == false);
        rhs = 100;
        expect(lhs >= rhs);
    };

    "=="_test = [] {
        gmprime lhs(100), rhs(50);
        expect((lhs == rhs) == false);
        rhs = 100;
        expect(lhs == rhs);
    };

    "!="_test = [] {
        gmprime lhs(100), rhs(50);
        expect(lhs != rhs);
        rhs = 100;
        expect((lhs != rhs) == false);
    };

    // Compared exactly, like an int is promoted to double, instead of truncating the floating point value
    "floating_point_fraction"_test = [] {
        const gmprime a(5), b(-5);

        expect(!(a == 5.5));
        expect(a != 5.5);
        expect(a < 5.5);
        expect(a <= 5.5);
        expect(!(a > 5.5));
        expect(!(a >= 5.5));
        expect(a > 4.5);
        expect(a == 5.0);
        expect(a <= 5.0 && a >= 5.0);

        expect(5.5 > a);
        expect(!(5.5 == a));
        expect(4.5 < a);

        expect(b < -4.5);
        expect(b > -5.5);
        expect(-5.5f < b);
    };
    "floating_point_limits"_test = [] {
        const gmprime zero, big("1" + std::string(400, '0'));
        const auto nan = std::numeric_limits<double>::quiet_NaN();
        const auto inf = std::numeric_limits<double>::infinity();

        expect(zero == -0.0);
        expect(zero < 1e-300);
        expect(zero > -1e-300);

        expect(big < inf);
        expect(-big > -inf);
        expect(inf > big);

        expect(gmprime(1e300) == 1e300) << "a double converted to gmprime compares equal to itself";
        expect(gmprime("1" + std::string(300, '0')) != 1e300) << "1e300 is not exactly 10^300";

        // NaN is unordered: every comparison is false, except !=
        expect(!(zero == nan));
        expect(zero != nan);
        expect(!(zero < nan) && !(zero <= nan) && !(zero > nan) && !(zero >= nan));
        expect(!(nan == zero) && nan != zero);
    };
    "floating_point_long_double"_test = [] {
        // Exact beyond the 53 bits of a double where long double is wider (e.g. x86 Linux, ARM64 Linux)
        if constexpr (std::numeric_limits<long double>::digits > 61) {
            const gmprime a = (gmprime(1) << 60) + 1;
            const long double x = std::ldexp(1.0L, 60) + 1;

            expect(a == x);
            expect(a - 1 < x);
        }
    };
};

suite<"gmprime::arithmetic"> arithmetic = [] {
    "+"_test = [] {
        const gmprime a(42), b(42);
        const gmprime r = a + b;
        expect(r == 84);
    };
    "+="_test = [] {
        gmprime a(42);
        const gmprime b(42);
        a += b;
        expect(a == 84);
    };
    "-"_test = [] {
        const gmprime a(42), b(42);
        const gmprime r = a - b;
        expect(r == 0);
    };
    "-="_test = [] {
        gmprime a(42);
        const gmprime b(42);
        a -= b;
        expect(a == 0);
    };
    "*"_test = [] {
        const gmprime a(42), b(42);
        const gmprime r = a * b;
        expect(r == 1764);
    };
    "*="_test = [] {
        gmprime a(42);
        const gmprime b(42);
        a *= b;
        expect(a == 1764);
    };
    "/"_test = [] {
        const gmprime a(42), b(42);
        const gmprime r = a / b;
        expect(r == 1);
    };
    "/="_test = [] {
        gmprime a(42);
        const gmprime b(42);
        a /= b;
        expect(a == 1);
    };
    "%"_test = [] {
        const gmprime a(42), b(42);
        const gmprime r = a % b;
        expect(r == 0);
    };
    "%="_test = [] {
        gmprime a(42);
        const gmprime b(42);
        a %= b;
        expect(a == 0);
    };
    "division_by_zero"_test = [] {
        gmprime a(42);
        const gmprime zero;

        expect(throws<std::domain_error>([&] { a / zero; }));
        expect(throws<std::domain_error>([&] { a % zero; }));
        expect(throws<std::domain_error>([&] { a /= zero; }));
        expect(throws<std::domain_error>([&] { a %= zero; }));
        expect(a == 42) << "operand modified by a failed division";
    };
    "++"_test = [] {
        gmprime r(42);
        ++r;
        expect(r == 43);
    };
    "--"_test = [] {
        gmprime r(42);
        --r;
        expect(r == 41);
    };
    "postfix_deleted"_test = [] {
        expect(!has_postfix_increment<gmprime>::value);
        expect(!has_postfix_decrement<gmprime>::value);
        expect(has_postfix_increment<int>::value) << "trait sanity check";
    };
    "unary+"_test = [] {
        const gmprime r(-42);
        expect(+r == -42);
    };
    "unary-"_test = [] {
        const gmprime r(42);
        expect(-r == -42);
    };

    "primitive_lhs"_test = [] {
        check_primitive_lhs(7);
        check_primitive_lhs(7L);
        check_primitive_lhs(7UL);
    };
};

suite<"gmprime::shift"> shift = [] {
    "<<"_test = [] {
        const gmprime a(42);
        const gmprime r = a << 3;
        expect(r == 336);
    };
    "<<="_test = [] {
        gmprime r(42);
        r <<= 3;
        expect(r == 336);
    };
    ">>"_test = [] {
        const gmprime a(336);
        const gmprime r = a >> 3;
        expect(r == 42);
    };
    ">>="_test = [] {
        gmprime r(336);
        r >>= 3;
        expect(r == 42);
    };
    "shift_negative"_test = [] {
        gmprime a(42);

        expect(throws<std::domain_error>([&] { static_cast<void>(a << -1); }));
        expect(throws<std::domain_error>([&] { static_cast<void>(a >> -1); }));
        expect(throws<std::domain_error>([&] { a <<= -1; }));
        expect(throws<std::domain_error>([&] { a >>= -1; }));
        expect(a == 42) << "operand modified by a failed shift";
    };
    "shift_left_beyond_capacity"_test = [] {
        gmprime a(42);

        // GMP aborts the program when a number exceeds its maximum size
        expect(throws<std::length_error>([&] { static_cast<void>(a << std::numeric_limits<size_t>::max()); }));
        expect(throws<std::length_error>([&] { a <<= std::numeric_limits<size_t>::max(); }));
        expect(a == 42) << "operand modified by a failed shift";
    };
    "shift_right_huge"_test = [] {
        // GMP takes shift counts as unsigned long, 32 bits on Windows: a huge count must not be truncated
        const auto huge = std::numeric_limits<size_t>::max();
        expect((gmprime(42) >> huge) == 0);
        expect((gmprime(-42) >> huge) == -1) << "right shift rounds toward minus infinity, like an int";

        gmprime a(42);
        a >>= huge;
        expect(a == 0);
    };
};

suite<"gmprime::logical"> logical = [] {
    "&"_test = [] {
        const gmprime a(42), b(24);
        const gmprime r = a & b;
        expect(r == 8);
    };
    "&="_test = [] {
        gmprime r(42);
        const gmprime b(24);
        r &= b;
        expect(r == 8);
    };
    "|"_test = [] {
        const gmprime a(42), b(24);
        const gmprime r = a | b;
        expect(r == 58);
    };
    "|="_test = [] {
        gmprime r(42);
        const gmprime b(24);
        r |= b;
        expect(r == 58);
    };
    "^"_test = [] {
        const gmprime a(42), b(24);
        const gmprime r = a ^ b;
        expect(r == 50);
    };
    "^="_test = [] {
        gmprime r(42);
        const gmprime b(24);
        r ^= b;
        expect(r == 50);
    };
    "~"_test = [] {
        const gmprime a(42);
        const gmprime r = ~a;
        constexpr int e = ~42;
        expect(r == e);
    };
    "explicit_bool"_test = [] {
        const gmprime positive(42), negative(-42), zero;

        expect(static_cast<bool>(positive));
        expect(static_cast<bool>(negative));
        expect(!zero);
        expect(positive && !zero);

        bool taken = false;
        if (positive) {
            taken = true;
        }
        expect(taken);

        expect(!std::is_convertible_v<gmprime, bool>) << "conversion to bool must be explicit";
    };
};

suite<"gmprime::accessor"> accessor = [] {
    "[]"_test = [] {
        gmprime r(ULONG_LONG_MAX);
        const mp_limb_t word = r[0];
        expect(word == ULONG_LONG_MAX);
    };
    "[]<unsigned char>"_test = [] {
        const gmprime r(4242);
        // const auto word0 = r.operator[]<unsigned char>(0);
        const auto word0 = r[0_ucw];
        // const auto word1 = r.operator[]<unsigned char>(1);
        const auto word1 = r[1_ucw];

        expect(word0 == 146);
        expect(word1 == 16);
    };
    "[]<unsigned char>"_test = [] {
        const gmprime r(4242);
        const auto word0 = r.operator[]<unsigned char>(0);
        // const unsigned char word0 = r[0];
        const auto word1 = r.operator[]<unsigned char>(1);
        // const unsigned char word1 = r[1];

        expect(word0 == 146);
        expect(word1 == 16);
    };
    "[]<unsigned char>="_test = [] {
        gmprime r(4242);
        const auto before = static_cast<unsigned char>(r.operator[]<unsigned char>(0));
        r.operator[]<unsigned char>(0) = 'a';
        const auto after = static_cast<unsigned char>(r.operator[]<unsigned char>(0));

        expect(before == 146);
        expect(after == 'a');
    };
    "[]<bool>"_test = [] {
        const gmprime r(4242);
        const bool bit1 = r.operator[]<bool>(0);
        const bool bit2 = r.operator[]<bool>(1);

        expect(bit1 == false);
        expect(bit2 == true);
    };
    "[]<bool>="_test = [] {
        gmprime r(4242);
        const auto before = static_cast<bool>(r.operator[]<bool>(0));
        // r.operator[]<bool>(0) = true;
        r[0_bw] = true;
        const auto after = static_cast<bool>(r.operator[]<bool>(0));

        expect(before == false);
        expect(after == true);
    };

    "[]<unsigned char>-="_test = [] {
        gmprime r(0x1000);
        r[1_ucw] -= 1;

        expect(r == 0xF00);
    };
    "[]<unsigned short>=_beyond_size"_test = [] {
        gmprime r(1);
        r.operator[]<unsigned short>(5) = 0xBEEF;

        expect(r == (gmprime(0xBEEF) << 80) + 1);
    };
    "[]<unsigned char>=_negative"_test = [] {
        gmprime r(-0x1092);
        r[0_ucw] = 0x61;

        expect(r == -0x1061) << "word write must keep the sign";
    };
    "[]<signed char>="_test = [] {
        gmprime r(0x1000);
        r.operator[]<signed char>(0) = -1;

        expect(r == 0x10FF);
        expect(r.operator[]<signed char>(0) == -1);
    };
#ifdef __SIZEOF_INT128__
    // A word wider than a limb spans several limbs
    "[]<unsigned __int128>"_test = [] {
        const gmprime r("0x0123456789abcdeffedcba98765432100011223344556677", 0);
        const auto word0 = r.operator[]<unsigned __int128>(0);
        const auto word1 = r.operator[]<unsigned __int128>(1);

        expect(gmprime(word0) == gmprime("0xfedcba98765432100011223344556677", 0));
        expect(gmprime(word1) == gmprime("0x0123456789abcdef", 0));
    };
    "[]<unsigned __int128>="_test = [] {
        gmprime r(1);
        r.operator[]<unsigned __int128>(1) = (static_cast<unsigned __int128>(0xAAAAAAAAAAAAAAAAULL) << 64) | 0x5555555555555555ULL;

        expect(r == gmprime("0xaaaaaaaaaaaaaaaa555555555555555500000000000000000000000000000001", 0));
    };
#endif
    "[]==_mixed_signedness"_test = [] {
        gmprime r(ULONG_LONG_MAX);

        expect(!(r[0] == -1)) << "an all-ones word is not equal to -1";
        expect(r[0] != -1);
        expect(r[0] == ULONG_LONG_MAX);

        gmprime s(42);
        expect(s[0] == 42);
        expect(s[0] != 41);
    };
    "[]<bool>_negative"_test = [] {
        // Bits follow the magnitude, like words do
        gmprime r(-3);

        expect(r[0_ucw] == 3);
        expect(r[0_bw] == true);
        expect(r[1_bw] == true);
        expect(r[2_bw] == false);
    };
    "set_bit_negative"_test = [] {
        gmprime r(-3);

        r.set_bit(2, true);
        expect(r == -7);

        r.set_bit(0, false);
        expect(r == -6);
    };
};

// Expected values computed independently (Python integers)
suite<"gmprime::bit_window"> bit_windows = [] {
    // 192 bits: three 64-bit limbs
    const gmprime n("0x0123456789abcdeffedcba98765432100011223344556677", 0);

    "read_unsigned"_test = [&] {
        const uint8_t byte = n.bits(4, 8);
        const uint32_t across_limbs = n.bits(48, 32);
        const uint64_t limb_sized = n.bits(32, 64);
        const uint16_t beyond_size = n.bits(1000, 13);

        expect(byte == 0x67);
        expect(across_limbs == 0x32100011u);
        expect(limb_sized == 0x7654321000112233ull);
        expect(beyond_size == 0);
    };
    "read_gmprime"_test = [&] {
        const gmprime wide = n.bits(40, 100);
        expect(wide == gmprime("0xdeffedcba9876543210001122", 0));
    };
    "read_negative_magnitude"_test = [&] {
        const gmprime m = -n;
        const uint8_t byte = m.bits(4, 8);
        const gmprime wide = m.bits(40, 100);

        expect(byte == 0x67);
        expect(wide == gmprime("0xdeffedcba9876543210001122", 0));
    };
    "read_wider_than_type"_test = [&] {
        expect(throws<std::overflow_error>([&] { [[maybe_unused]] const uint8_t x = n.bits(0, 9); }));
        expect(throws<std::overflow_error>([&] { [[maybe_unused]] const uint64_t x = n.bits(0, 65); }));
    };
    "read_empty"_test = [&] {
        const uint8_t x = n.bits(10, 0);
        const gmprime g = n.bits(10, 0);
        expect(x == 0);
        expect(g == 0);
    };
    "word_access_is_a_window"_test = [&] {
        const uint16_t window = n.bits(5 * 16, 16);
        expect(n.operator[]<uint16_t>(5) == window);
        expect(window == 0x7654);
    };
    "ostream"_test = [&] {
        std::stringstream t;
        t << n.bits(4, 8) << ' ' << gmprime::hex << n.bits(40, 100) << ' ' << gmprime::bin << n.bits(4, 8);
        expect(t.str() == "103 deffedcba9876543210001122 1100111");
    };
    "position_width"_test = [&] {
        const auto w = n.bits(1000, 13);
        expect(w.position() == 1000);
        expect(w.width() == 13);
    };

    "write_unsigned"_test = [&] {
        gmprime a(n);
        a.bits(4, 8) = 0xABu;
        expect(a == gmprime("0x123456789abcdeffedcba98765432100011223344556ab7", 0));
    };
    "write_across_limbs"_test = [&] {
        gmprime a(n);
        a.bits(60, 8) = 0xFFu;
        expect(a == gmprime("0x123456789abcdeffedcba987654321ff011223344556677", 0));
    };
    "write_truncates"_test = [] {
        gmprime a;
        a.bits(0, 4) = 0xFFu;
        expect(a == 0xF);

        gmprime b;
        b.bits(8, 4) = gmprime(0x1FF);
        expect(b == 0xF00);
    };
    "write_beyond_size"_test = [] {
        gmprime a(1);
        a.bits(100, 8) = 0x5Au;
        expect(a == (gmprime(0x5A) << 100) + 1);
    };
    "write_keeps_sign"_test = [&] {
        gmprime a = -n;
        a.bits(4, 8) = 0xABu;
        expect(a == gmprime("-0x123456789abcdeffedcba98765432100011223344556ab7", 0));

        // Wider than any integer: goes through GMP operations instead of limb masks
        gmprime b = -n;
        b.bits(40, 100) = gmprime("0x0f0f0f0f0f0f0f0f0f0f0f0f0f", 0);
        expect(b == gmprime("-0x123456789abcf0f0f0f0f0f0f0f0f0f0f0f0f3344556677", 0));
    };
    "write_gmprime"_test = [&] {
        gmprime a(n);
        a.bits(40, 100) = gmprime("0x0f0f0f0f0f0f0f0f0f0f0f0f0f", 0);
        expect(a == gmprime("0x123456789abcf0f0f0f0f0f0f0f0f0f0f0f0f3344556677", 0));

        gmprime b;
        b.bits(0, 8) = gmprime(-5);
        expect(b == 5) << "the magnitude of the value is written";
    };
    "write_zero_extends"_test = [] {
        gmprime a = (gmprime(1) << 200) - 1;
        a.bits(0, 100) = static_cast<uint8_t>(5);
        expect(a == gmprime("0xfffffffffffffffffffffffff0000000000000000000000005", 0));
    };
    "write_empty"_test = [&] {
        gmprime a(n);
        a.bits(10, 0) = 0xFFu;
        expect(a == n);
    };
    "write_window_copies_bits"_test = [] {
        gmprime a(0x1234);
        a.bits(0, 8) = a.bits(8, 8);
        expect(a == 0x1212) << "assigning a window copies its bits, it does not rebind";

        const gmprime c(0xAB);
        a.bits(8, 8) = c.bits(0, 8);
        expect(a == 0xAB12);
    };

    "unsigned_only"_test = [] {
        expect(std::is_assignable_v<bit_window<false>, unsigned>);
        expect(std::is_assignable_v<bit_window<false>, const gmprime&>);
        expect(!std::is_assignable_v<bit_window<false>, int>) << "signed integers are refused";
        expect(!std::is_assignable_v<bit_window<false>, bool>);
        expect(!std::is_assignable_v<bit_window<false>, char>);
        expect(!std::is_assignable_v<bit_window<false>, double>) << "would silently drop the sign through gmprime";

        expect(std::is_convertible_v<bit_window<false>, unsigned>);
        expect(std::is_convertible_v<bit_window<false>, gmprime>);
        expect(!std::is_convertible_v<bit_window<false>, int>) << "signed integers are refused";
        expect(!std::is_convertible_v<bit_window<false>, char>);
    };
    "const_is_read_only"_test = [] {
        expect(std::is_convertible_v<bit_window<true>, unsigned>);
        expect(std::is_convertible_v<bit_window<true>, gmprime>);
        expect(!std::is_assignable_v<bit_window<true>, unsigned>);
        expect(!std::is_assignable_v<bit_window<true>, const gmprime&>);
    };
};

// Each function is also checked against a naive bit-by-bit version on large irregular numbers
suite<"gmprime::bit_analysis"> bit_analysis = [] {
    // Bits 3, 64 and 200: in the first, second and fourth limbs
    const gmprime sparse = (gmprime(1) << 3) + (gmprime(1) << 64) + (gmprime(1) << 200);
    const gmprime factorial = gmprime::factorial(300);
    const std::vector<gmprime> irregular = {
        factorial,
        factorial.next_prime(),
        (gmprime(1) << 1000) - 1,
        gmprime("0x0123456789abcdeffedcba9876543210", 0),
    };

    "popcount"_test = [&] {
        expect(gmprime().popcount() == 0);
        expect(gmprime(0xFF).popcount() == 8);
        expect(sparse.popcount() == 3);
        expect(((gmprime(1) << 1000) - 1).popcount() == 1000);
        expect(gmprime(-0xFF).popcount() == 8) << "the magnitude, not an infinity of two's complement bits";

        for (const auto& value : irregular) {
            std::vector<size_t> expected;
            naive_set_bits(value, expected);
            expect(value.popcount() == expected.size());
        }
    };
    "popcount_text_file"_test = [] {
        const auto path = TEST_OUTPUT_DIR "popcount.txt";
        {
            std::ofstream out(path, std::ios::binary);
            out << "abc";
        }

        gmprime text;
        text.load_from_file(path);
        expect(text.popcount() == 10) << "a = 0x61, b = 0x62, c = 0x63: 3 + 3 + 4 bits set";

        std::remove(path);
    };

    "hamming_distance"_test = [&] {
        expect(gmprime(0b1010).hamming_distance(gmprime(0b0110)) == 2);
        expect(gmprime(1).hamming_distance(sparse) == 4) << "operands of different sizes: bits 0, 3, 64 and 200";
        expect(sparse.hamming_distance(gmprime(1)) == 4) << "symmetric";
        expect(gmprime().hamming_distance(sparse) == sparse.popcount());
        expect(gmprime(-5).hamming_distance(gmprime(5)) == 0) << "on the magnitudes";

        for (const auto& a : irregular) {
            for (const auto& b : irregular) {
                expect(a.hamming_distance(b) == (a ^ b).popcount());
            }
        }
    };

    "next_set_bit"_test = [&] {
        expect(sparse.next_set_bit() == 3);
        expect(sparse.next_set_bit(3) == 3);
        expect(sparse.next_set_bit(4) == 64) << "skips the rest of a limb";
        expect(sparse.next_set_bit(65) == 200) << "skips whole zero limbs";
        expect(sparse.next_set_bit(201) == gmprime::npos);
        expect(sparse.next_set_bit(100000) == gmprime::npos);
        expect(gmprime().next_set_bit() == gmprime::npos);
        expect(gmprime(-8).next_set_bit() == 3) << "on the magnitude";
    };

    "next_clear_bit"_test = [] {
        expect(gmprime(0b0111).next_clear_bit() == 3);
        expect(gmprime(0b0111).next_clear_bit(5) == 5);
        expect(((gmprime(1) << 64) - 1).next_clear_bit() == 64) << "beyond the magnitude, bits are clear";
        expect(((gmprime(1) << 130) - 1).next_clear_bit(10) == 130) << "skips whole full limbs";
        expect(((gmprime(1) << 130) - 1).next_clear_bit(1000) == 1000);
        expect(gmprime().next_clear_bit(7) == 7);
        expect(gmprime(-7).next_clear_bit() == 3) << "on the magnitude";
    };

    "set_bits"_test = [&] {
        const std::vector<size_t> positions(sparse.set_bits().begin(), sparse.set_bits().end());
        expect(positions == std::vector<size_t>{3, 64, 200});

        const auto empty = gmprime().set_bits();
        expect(empty.begin() == empty.end());

        for (const auto& value : irregular) {
            std::vector<size_t> expected;
            naive_set_bits(value, expected);
            const auto range = value.set_bits();
            expect(std::vector<size_t>(range.begin(), range.end()) == expected);
        }
    };
    "set_bits_temporary"_test = [&] {
        // The temporary would be destroyed before the loop runs if the range only pointed to it
        std::vector<size_t> positions;
        for (const size_t pos : (sparse | gmprime(1)).set_bits()) {
            positions.push_back(pos);
        }
        expect(positions == std::vector<size_t>{0, 3, 64, 200});

        // Copies and moves of an owning range keep a valid number
        auto owning = (sparse & gmprime(0xFF)).set_bits();
        const auto copy = owning;
        const auto moved = std::move(owning);
        expect(std::vector<size_t>(copy.begin(), copy.end()) == std::vector<size_t>{3});
        expect(std::vector<size_t>(moved.begin(), moved.end()) == std::vector<size_t>{3});
    };
    "set_bits_iterator"_test = [&] {
        const auto range = sparse.set_bits();
        auto it = range.begin();

        expect(*it == 3);
        expect(*it++ == 3) << "postfix returns the previous position";
        expect(*it == 64);
        expect(std::distance(range.begin(), range.end()) == 3);
        expect(std::is_same_v<std::iterator_traits<set_bit_iterator>::iterator_category, std::forward_iterator_tag>);
    };

    "histogram_small"_test = [] {
        // 0b1101 over 4 bits, sliding windows of 2 bits: bits 0-1 hold 01 (value 1), bits 1-2 hold 10
        // (value 2, bit 1 being its low bit), bits 2-3 hold 11 (value 3)
        expect(gmprime(0b1101).histogram(2, 1, 4) == std::vector<size_t>{0, 1, 1, 1});
        expect(gmprime(0b1101).histogram(2, 0, 4) == std::vector<size_t>{0, 1, 0, 1}) << "aligned: bits 0-1, then 2-3";
        expect(gmprime(0x0F0F).histogram(4, 8) == std::vector<size_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2})
            << "step wider than width: the low nibble of each byte";
        expect(gmprime(1).histogram(8, 0, 24)[0] == 2) << "explicit length: the zero bytes above count";
        expect(gmprime(-0b1101).histogram(2) == gmprime(0b1101).histogram(2)) << "on the magnitude";
        expect(gmprime(0b1101).histogram(8, 0, 4) == std::vector<size_t>(256, 0)) << "no complete window";
        expect(gmprime(0b1101).histogram(8)[0b1101] == 1) << "by default, the number in whole bytes, like a file";
        expect(gmprime(0b1101).histogram(2)[0] == 2) << "so the high bits of its partial byte count";
    };
    "histogram_invalid_width"_test = [] {
        expect(throws<std::invalid_argument>([] { static_cast<void>(gmprime(1).histogram(0)); }));
        constexpr auto too_wide = gmprime::max_histogram_width + 1;
        expect(throws<std::invalid_argument>([] { static_cast<void>(gmprime(1).histogram(too_wide)); }));
    };
    "histogram_matches_windows"_test = [&] {
        for (const auto& value : irregular) {
            const size_t length = (value.used_bits() + 7) / 8 * 8;
            for (const size_t width : {1, 3, 8, 13}) {
                for (const size_t step : {size_t{1}, width, size_t{5}}) {
                    std::vector<size_t> expected(size_t{1} << width, 0);
                    for (size_t pos = 0; pos + width <= length; pos += step) {
                        const uint32_t symbol = value.bits(pos, width);
                        ++expected[symbol];
                    }
                    expect(value.histogram(width, step) == expected);
                }
            }
        }
    };
    "histogram_text_file"_test = [] {
        const auto path = TEST_OUTPUT_DIR "histogram.txt";
        gmprime text;

        {
            std::ofstream out(path, std::ios::binary);
            out << "abca";
        }
        text.load_from_file(path);
        const auto bytes = text.histogram(8);
        expect(bytes['a'] == 2 && bytes['b'] == 1 && bytes['c'] == 1);

        {
            std::ofstream out(path, std::ios::binary);
            out.write("ab\0\0", 4);
        }
        text.load_from_file(path);
        expect(text.histogram(8)[0] == 2) << "the file's high zero bytes count, though GMP does not store them";

        std::remove(path);
    };

    "entropy"_test = [] {
        const auto close = [](const double a, const double b) { return std::fabs(a - b) < 1e-12; };

        expect(close(gmprime::entropy({}), 0));
        expect(close(gmprime::entropy({0, 0}), 0));
        expect(close(gmprime::entropy({5}), 0)) << "a single symbol carries no information";
        expect(close(gmprime::entropy({7, 7}), 1));
        expect(close(gmprime::entropy({3, 1}), 0.8112781244591328));
        expect(close(gmprime::entropy(std::vector<size_t>(256, 1)), 8)) << "uniform bytes: 8 bits per byte";
    };
    "entropy_text"_test = [] {
        const gmprime same("0x61616161", 0);
        expect(gmprime::entropy(same.histogram(8)) == 0) << "aaaa";

        // Every byte value once: incompressible byte by byte
        gmprime all;
        for (unsigned value = 0; value < 256; ++value) {
            all.bits(8 * value, 8) = value;
        }
        expect(std::fabs(gmprime::entropy(all.histogram(8)) - 8) < 1e-12);
    };
};

// Expected values computed independently (Python), or checked against histogram() and naive versions
suite<"gmprime::context_analysis"> context_analysis = [] {
    const auto close = [](const double a, const double b) { return std::fabs(a - b) < 1e-9; };

    "conditional_entropy_order_0"_test = [&] {
        std::string text;
        random_text(500, 7, 42, text);
        gmprime value;
        from_bytes(text, value);

        for (const size_t width : {1, 3, 8, 13}) {
            expect(close(value.conditional_entropy(0, width).naive, gmprime::entropy(value.histogram(width))));
        }
    };
    "conditional_entropy_matches_histograms"_test = [&] {
        std::string text;
        random_text(2000, 5, 7, text);
        gmprime value;
        from_bytes(text, value);
        const size_t n = text.size();

        // H(X | previous symbols) = H(windows of order + 1 symbols) - H(windows of order symbols), over the
        // same starts
        expect(close(value.conditional_entropy(1).naive,
                     gmprime::entropy(value.histogram(16, 8)) - gmprime::entropy(value.histogram(8, 8, 8 * (n - 1)))));
        expect(close(value.conditional_entropy(2).naive,
                     gmprime::entropy(value.histogram(24, 8)) - gmprime::entropy(value.histogram(16, 8, 8 * (n - 1)))));
    };
    "conditional_entropy_values"_test = [&] {
        gmprime value;
        from_bytes("abababab", value);
        const auto alternating = value.conditional_entropy(1);
        expect(alternating.naive == 0 && alternating.corrected == 0) << "each symbol determined by the previous one";
        expect(alternating.contexts == 2 && alternating.samples == 7);

        gmprime small;
        from_bytes("aab", small);
        const auto skewed = small.conditional_entropy(0);
        expect(close(skewed.naive, 0.9182958340544896));
        expect(close(skewed.corrected, 1.1587450075359835)) << "Miller-Madow: + (2 - 1) / (2 * 3 * ln 2)";

        expect(gmprime().conditional_entropy(2).samples == 0) << "fewer symbols than the order";
        expect(throws<std::invalid_argument>([&] { static_cast<void>(value.conditional_entropy(1, 0)); }));
        expect(throws<std::invalid_argument>([&] { static_cast<void>(value.conditional_entropy(1, 33)); }));
    };

    "adaptive_code_length_values"_test = [&] {
        gmprime ab, aaaa;
        from_bytes("ab", ab);
        from_bytes("aaaa", aaaa);

        expect(close(ab.adaptive_code_length(0), 16.011227255423254));
        expect(close(aaaa.adaptive_code_length(1), 28.126704472843187));
        expect(gmprime().adaptive_code_length(0) == 0);
        expect(ab.adaptive_code_length(5) == 16) << "every symbol in the uncoded prefix: width bits each";
    };
    "adaptive_code_length_above_entropy"_test = [&] {
        std::string text;
        random_text(3000, 6, 3, text);
        gmprime value;
        from_bytes(text, value);

        // The empirical entropy is the shortest code length knowing the frequencies in advance
        for (const size_t order : {0, 1, 2}) {
            const auto entropy = value.conditional_entropy(order);
            const auto bits = value.adaptive_code_length(order);
            expect(bits > entropy.naive * static_cast<double>(entropy.samples));
        }
    };

    "longest_previous_factors_values"_test = [] {
        const auto lpf = [](const std::string& text) {
            gmprime value;
            from_bytes(text, value);
            return value.longest_previous_factors();
        };

        expect(lpf("abcabc") == std::vector<uint32_t>{0, 0, 0, 3, 2, 1});
        expect(lpf("aaaa") == std::vector<uint32_t>{0, 3, 2, 1}) << "a factor may overlap its occurrence";
        expect(lpf("abab") == std::vector<uint32_t>{0, 0, 2, 1});
        expect(lpf("a") == std::vector<uint32_t>{0});
        expect(gmprime().longest_previous_factors().empty());
    };
    "longest_previous_factors_file"_test = [] {
        const auto path = TEST_OUTPUT_DIR "repetitions.bin";
        {
            std::ofstream out(path, std::ios::binary);
            out.write("ab\0\0ab\0\0", 8);
        }
        gmprime value;
        value.load_from_file(path);

        expect(value.longest_previous_factors() == std::vector<uint32_t>{0, 0, 0, 1, 4, 3, 2, 1})
            << "the file's high zero bytes count";

        std::remove(path);
    };
    "longest_previous_factors_matches_naive"_test = [] {
        for (const unsigned alphabet : {2u, 3u, 26u}) {
            for (const size_t length : {2, 3, 5, 8, 13, 40, 100}) {
                for (uint32_t seed = 1; seed <= 20; ++seed) {
                    std::string text;
                    random_text(length, alphabet, seed, text);
                    gmprime value;
                    from_bytes(text, value);
                    std::vector<uint32_t> expected;
                    naive_longest_previous_factors(text, expected);

                    expect(value.longest_previous_factors() == expected) << text;
                }
            }
        }
    };

    "repetitions"_test = [] {
        gmprime abc, aaaa;
        from_bytes("abcabcabc", abc);
        from_bytes("aaaa", aaaa);

        const auto three = abc.repetitions(3);
        expect(three.bytes == 9 && three.matches == 1 && three.matched_bytes == 6 && three.literals == 3);

        const auto none = abc.repetitions(7);
        expect(none.matches == 0 && none.literals == 9) << "no repeat that long";

        const auto any = aaaa.repetitions(0);
        expect(any.matches == 1 && any.matched_bytes == 3 && any.literals == 1) << "min_length 0 counts as 1";
    };

    "bytes_round_trip"_test = [] {
        const std::vector<unsigned char> data{'a', 0, 'b', 0, 0};
        const auto value = gmprime::from_bytes(data);

        expect(value.bytes() == data) << "high zero bytes kept, like a loaded file";
        expect(value == 0x620061);
        expect(gmprime::from_bytes({}).bytes().empty());
    };

    "autocorrelation"_test = [] {
        // "ab" repeated: the bits repeat every 16
        const auto periodic = gmprime::from_bytes(std::vector<unsigned char>{'a', 'b', 'a', 'b', 'a', 'b', 'a', 'b'});
        const auto self = periodic.autocorrelation(100);

        expect(self.size() == 101);
        expect(self[0] == 64) << "every bit equals itself";
        expect(self[16] == 48) << "a period: all of the overlap equal";
        expect(self[8] < 48);
        expect(self[64] == 0 && self[100] == 0) << "no overlap";
    };
    "autocorrelation_matches_naive"_test = [] {
        std::string text;
        random_text(200, 4, 11, text);
        gmprime value;
        from_bytes(text, value);
        const size_t length = 8 * text.size();
        const auto self = value.autocorrelation(70);

        for (size_t d = 0; d <= 70; ++d) {
            size_t equal = 0;
            for (size_t p = 0; p + d < length; ++p) {
                equal += value.get_bit(p) == value.get_bit(p + d) ? 1 : 0;
            }
            expect(self[d] == equal) << "distance" << d;
        }
        expect(gmprime(-0x1234).autocorrelation(5) == gmprime(0x1234).autocorrelation(5)) << "on the magnitude";
    };

    "stream_order"_test = [] {
        const std::vector<unsigned char> data{0b00000001, 0b11000000, 0};
        const auto stream = gmprime::from_bytes(data).stream_order();

        expect(stream.bytes() == std::vector<unsigned char>{0b10000000, 0b00000011, 0}) << "each byte reversed, size kept";
        for (size_t p = 0; p < 24; ++p) {
            expect(static_cast<int>(stream.get_bit(p)) == GMPrime::coding::stream_bit(data, p)) << "bit" << p;
        }
    };

    "analyze_folds"_test = [] {
        // A bitmap of 60 rows of 200 bits (not a multiple of 8), each row the previous one with a few flips
        const size_t width = 200, rows = 60;
        std::vector<unsigned char> bitmap(width * rows / 8, 0);
        uint32_t seed = 7;
        std::vector<int> row(width);
        for (auto& bit : row) {
            seed = seed * 1664525u + 1013904223u;
            bit = (seed >> 20) & 1;
        }
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < width; ++c) {
                seed = seed * 1664525u + 1013904223u;
                if ((seed >> 16) % 64 == 0) row[c] ^= 1;
                const size_t p = r * width + c;
                bitmap[p / 8] = static_cast<unsigned char>(bitmap[p / 8] | (row[c] << (7 - p % 8)));
            }
        }
        const auto folds = gmprime::from_bytes(bitmap).analyze_folds();
        expect(folds.width == width) << "found the row width:" << folds.width;
        expect(folds.strong());
        expect(folds.score < folds.median && (folds.layer == 0 || folds.layer >= 3));

        std::string text;
        random_text(3000, 26, 9, text);
        gmprime value;
        from_bytes(text, value);
        expect(!value.analyze_folds().strong()) << "text without rows";

        const auto tiny = gmprime(5).analyze_folds();
        expect(tiny.width == 0 && !tiny.strong()) << "too short to search";
        expect(value.analyze_folds(100, 50).width == 0) << "no width to search";

        // A sample of 40 rows finds the same fold; a sample larger than the data analyzes all of it
        const auto sampled = gmprime::from_bytes(bitmap).analyze_folds(64, 8192, width * 40 / 8);
        expect(sampled.width == width && sampled.strong());
        const auto all = gmprime::from_bytes(bitmap).analyze_folds(64, 8192, bitmap.size() * 2);
        expect(all.width == folds.width && all.score == folds.score && all.layer == folds.layer);
    };
};

suite<"gmprime::coding"> coding_suite = [] {
    using namespace GMPrime::coding;

    "logistic_domain"_test = [] {
        expect(squash(0) == 32768);
        expect(squash(3072) == 65535 && squash(5000) == 65535 && squash(-3072) == 1 && squash(-5000) == 1) << "clamped beyond 12";

        // Within what the 12-bit coder can tell apart of 65536 / (1 + e^(-d / 256))
        double worst_squash = 0;
        int worst_round_trip = 0;
        for (int d = -3071; d <= 3071; ++d) {
            worst_squash = std::max(worst_squash, std::fabs(squash(d) - 65536 / (1 + std::exp(-d / 256.0))));
            if (d >= -1500 && d <= 1500) {
                worst_round_trip = std::max(worst_round_trip, std::abs(stretch(squash(d)) - d));
            }
        }
        expect(worst_squash < 16) << "squash error" << worst_squash;
        expect(worst_round_trip <= 1) << "stretch inverts squash" << worst_round_trip;
        expect(stretch(32768) == 0 && stretch(-5) == stretch(0) && stretch(70000) == stretch(65535));

        expect(clamp_probability(0) == 16 && clamp_probability(65535) == 65520 && clamp_probability(1000) == 1000);
        expect(probability_12(32768) == 2048 && probability_12(0) == 1 && probability_12(65535) == 4095);
    };

    "counter_rate"_test = [] {
        // 65536 / (n + 1.5), by integer division
        expect(counter_rate(0) == 43690 && counter_rate(1) == 26214 && counter_rate(1023) == 63);
    };

    "bit_counter"_test = [] {
        bit_counter counter;
        expect(!counter.seen() && counter.p16() == 32768);

        counter.update(1, 3);
        expect(counter.seen() && counter.p16() > 52428) << "learns fast at first: rate 1 / 1.5";
        for (int i = 0; i < 10; ++i) {
            counter.update(1, 3);
        }
        expect(counter.n == 3) << "the count stops at the limit";
        expect(counter.p16() > 62259);
    };

    "bit_history"_test = [] {
        bit_history history;
        for (int i = 0; i < 3; ++i) {
            history.update(1);
        }
        expect(history.ones() == 3 && history.zeros() == 0);

        history.update(0);
        expect(history.zeros() == 1 && history.ones() == 2) << "the opposite count beyond 2 is halved, rounded up";

        for (int i = 0; i < 20; ++i) {
            history.update(0);
        }
        expect(history.zeros() == 15 && history.ones() == 2) << "capped at 15, the opposite count no longer halved at 2";
    };

    "history_map"_test = [] {
        history_map map;
        bit_history seen_one;
        seen_one.update(1);
        expect(std::abs(map.p(seen_one) - 49152) <= 1) << "starts from the counts: (1 + 1/2) / (0 + 1 + 1)";

        for (int i = 0; i < 50; ++i) {
            map.update(seen_one, 0);
        }
        expect(map.p(seen_one) < 13107) << "learns what follows that history in practice";
    };

    "hashed_tables"_test = [] {
        hashed_counters counters(10);
        counters.at(42).update(1, 30);
        expect(counters.at(42).seen()) << "same hash, same slot";

        const std::vector<uint8_t> history{'a', 'b', 'c'};
        expect(context_hash(history, 3, 1) != context_hash(history, 3, 2));
        expect(context_hash(history, 3, 0) != context_hash(history, 2, 1));
        expect(context_hash(history, 1, 2) == context_hash(std::vector<uint8_t>{0, 'a'}, 2, 2)) << "missing bytes count as 0";
    };

    "checked_table"_test = [] {
        // A 2-slot table: find two hashes sharing a slot with different checksums
        checked_table<context_slot> table(1);
        uint64_t other = 0;
        do {
            ++other;
        } while ((hash_mix(other) & 1) != (hash_mix(0) & 1) || ((hash_mix(other) >> 56) | 1) == ((hash_mix(0) >> 56) | 1));

        table.at(0).update(1, 30);
        expect(table.at(0).seen()) << "same context, same statistics";

        expect(!table.at(other).seen()) << "another context resets the slot instead of sharing it";
        expect(!table.at(0).seen()) << "and the newest context owns it";
    };
    "checked_table_buckets"_test = [] {
        checked_table<context_slot> table(10);
        table.prefetch(42);
        expect(!table.at(42).seen()) << "a prefetch changes nothing";

        // The 16 offsets of one hash are distinct slots, all owned by that hash
        for (size_t offset = 0; offset < 16; ++offset) {
            table.at(7, offset).update(1, 30);
            table.at(7, offset).n = static_cast<uint16_t>(offset + 1);
        }
        bool distinct = true;
        for (size_t offset = 0; offset < 16; ++offset) {
            distinct = distinct && table.at(7, offset).n == offset + 1;
        }
        expect(distinct);
    };

    "context_slot"_test = [] {
        context_slot slot;
        expect(!slot.seen());
        slot.update(1, 30);
        slot.update(1, 30);
        expect(slot.seen() && slot.n == 2 && slot.p16() > 58900) << "1/2, then 5/6, then 9/10";
        expect(slot.history.ones() == 2 && slot.history.zeros() == 0);
        expect(sizeof(context_slot) == 8);
    };

    "last_byte_slot"_test = [] {
        last_byte_slot slot;
        expect(slot.expected_bit(1, 0) == -1) << "nothing remembered";

        slot.update('a');
        slot.update('a');
        expect(slot.byte == 'a' && slot.count == 2);
        // 'a' = 0b01100001: first bit 0, then 1 once 0 was coded
        expect(slot.expected_bit(1, 0) == 0);
        expect(slot.expected_bit(0b10, 1) == 1);
        expect(slot.expected_bit(0b11, 1) == -1) << "the byte coded so far differs from 'a'";

        slot.update('b');
        expect(slot.byte == 'b' && slot.count == 1);
        for (int i = 0; i < 300; ++i) {
            slot.update('b');
        }
        expect(slot.count == 255);
    };

    "logistic_mixer"_test = [] {
        // Inputs in 1/256 units: 256 is 1, 77 about 0.3; rate 3277 is 0.05
        logistic_mixer mixer(2, 2, 3277);
        for (int i = 0; i < 200; ++i) {
            mixer.add(256);
            mixer.add(77);
            mixer.mix(0);
            mixer.update(1);
        }
        mixer.add(256);
        mixer.add(77);
        expect(mixer.mix(0) > 58982) << "learns to trust an input that predicts well";
        mixer.update(1);

        mixer.add(256);
        mixer.add(77);
        // Initial weights 19661 / 65536, about 0.3, in 8.24
        expect(mixer.mix(1) == squash(static_cast<int>((int64_t{19661} * 256 * (256 + 77)) >> 24))) << "other sets keep their initial weights";
        mixer.update(1);

        mixer.set_rate(0);
        mixer.add(256);
        mixer.add(77);
        const int before = mixer.mix(1);
        mixer.update(0);
        mixer.add(256);
        mixer.add(77);
        expect(mixer.mix(1) == before) << "a zero rate freezes the weights";
        mixer.update(0);

        logistic_mixer saturated(1, 1, 3277, 65536 * 4);
        saturated.add(3071);
        const std::vector<int> negative{-3071};
        expect(saturated.mix(0) == 65520 && saturated.mix(negative, 0) == 16) << "clamped for the coder";
        saturated.update(0);
    };

    "logistic_mixer_shared_inputs"_test = [] {
        // The same inputs given in place or added one by one mix and learn the same way
        logistic_mixer shared(3, 2, 3277), added(3, 2, 3277);
        const std::vector<int> inputs{128, -256, 512};
        for (int i = 0; i < 30; ++i) {
            for (const int x : inputs) {
                added.add(x);
            }
            const int p = shared.mix(inputs, i % 2);
            expect(p == added.mix(i % 2));
            shared.update(i % 3 == 0);
            added.update(i % 3 == 0);
        }

        // A copy of a mixer holding its own inputs learns from its copy of them
        logistic_mixer original(2, 1, 3277);
        original.add(256);
        original.add(128);
        const int before = original.mix(0);
        logistic_mixer copy = original;
        original.update(1);
        copy.update(1);
        original.add(256);
        original.add(128);
        copy.add(256);
        copy.add(128);
        expect(original.mix(0) == copy.mix(0) && copy.mix(0) > before);
    };

    "probability_map"_test = [] {
        // 45875 is about 0.7; rate 6554 is 0.1
        probability_map map(2, 6554);
        expect(std::abs(map.refine(45875, 0) - 45875) < 1311) << "starts close to the identity";

        for (int i = 0; i < 200; ++i) {
            map.refine(45875, 0);
            map.update(0);
        }
        expect(map.refine(45875, 0) < 13107) << "learns that 0.7 means 0 in context 0";
        map.update(0);
        expect(std::abs(map.refine(45875, 1) - 45875) < 1311) << "context 1 unchanged";
        map.update(1);
    };

    "match_model"_test = [] {
        match_model match(3, 10);
        std::vector<uint8_t> history;
        for (const char c : std::string("abcdef abcd")) {
            history.push_back(static_cast<uint8_t>(c));
            match.byte_added(history);
        }

        expect(match.length() == 4) << "abcd seen before, matched back to its start";
        // Predicts 'e' = 0b01100101 next
        expect(match.expected_bit(history, 1, 0) == 0);
        expect(match.expected_bit(history, 0b1011, 3) == 0) << "011 coded so far: e's prefix, next bit 0";
        expect(match.expected_bit(history, 0b1010, 3) == -1) << "010 coded so far: not e's prefix";

        match.miss();
        expect(match.length() == 0 && match.expected_bit(history, 1, 0) == -1);
    };

    "arithmetic_coder"_test = [] {
        std::string text;
        random_text(2000, 3, 5, text);

        // An adaptive order-1 bit model, rebuilt identically by the decoder
        const auto code = [&](auto&& code_bit) {
            hashed_counters counters(16);
            uint8_t previous = 0;
            std::vector<uint8_t> out;
            for (size_t n = 0; n < text.size(); ++n) {
                uint32_t c0 = 1;
                for (int i = 7; i >= 0; --i) {
                    auto& counter = counters.at(previous * 256u + c0);
                    const int bit = code_bit(clamp_probability(counter.p16()), (static_cast<uint8_t>(text[n]) >> i) & 1);
                    counter.update(bit, 255);
                    c0 = (c0 << 1) | static_cast<uint32_t>(bit);
                }
                previous = static_cast<uint8_t>(c0 & 0xFF);
                out.push_back(previous);
            }
            return out;
        };

        arithmetic_encoder encoder;
        code([&](const int p, const int bit) { encoder.encode(bit, p); return bit; });
        const gmprime compressed = encoder.finish();

        arithmetic_decoder decoder(compressed);
        const auto decoded = code([&](const int p, int) { return decoder.decode(p); });

        expect(decoded == std::vector<uint8_t>(text.begin(), text.end())) << "decodes what was encoded";
        expect(compressed.bytes().size() < text.size() / 2) << "3 letters: under 2 bits each";
    };
    "arithmetic_coder_extremes"_test = [] {
        // Certain probabilities, right and wrong, and an empty message
        arithmetic_encoder encoder;
        const int bits[] = {1, 1, 0, 1, 0, 0, 1};
        const int ps[] = {65470, 66, 32768, 65535, 0, 65530, 32768};
        for (int i = 0; i < 7; ++i) {
            encoder.encode(bits[i], ps[i]);
        }
        const gmprime compressed = encoder.finish();

        arithmetic_decoder decoder(compressed);
        bool same = true;
        for (int i = 0; i < 7; ++i) {
            same = same && decoder.decode(ps[i]) == bits[i];
        }
        expect(same);

        arithmetic_encoder empty;
        expect(empty.finish().bytes().size() == 4) << "the flush only";
    };

    "classify"_test = [] {
        expect(classify('\n') == char_class::newline && classify('\r') == char_class::newline);
        expect(classify(' ') == char_class::space && classify('\t') == char_class::space);
        expect(classify('q') == char_class::lowercase && classify('Q') == char_class::uppercase);
        expect(classify('7') == char_class::digit && classify(',') == char_class::punctuation);
        expect(classify(0) == char_class::other && classify(0xE9) == char_class::other);
    };

    "line_tracker"_test = [] {
        line_tracker lines;
        std::vector<uint8_t> history;
        const auto add = [&](const std::string& text) {
            for (const char c : text) {
                history.push_back(static_cast<uint8_t>(c));
                lines.byte_added(history);
            }
        };

        add("abc");
        expect(lines.column(history) == 3 && lines.above(history) == 0) << "no line above yet";

        add("defg\nxy");
        expect(lines.column(history) == 2);
        expect(lines.previous_length() == 8) << "abcdefg and its newline";
        expect(lines.above(history) == 'c') << "same column, previous line";

        add("zzzzzzzzzz");
        expect(lines.above(history) == 0) << "beyond the previous line";
    };

    "capitals"_test = [] {
        const std::string text = "The END. I saw NASA, McDonald and iPhone; the end\n";
        const std::vector<uint8_t> original(text.begin(), text.end());
        const auto encoded = capitals_encode(original);

        expect(encoded.has_value());
        const std::string lowered(encoded->begin(), encoded->end());
        expect(lowered == "\x01the \x02""end. \x01i saw \x02nasa, McDonald and iPhone; the end\n") << lowered;
        expect(capitals_decode(*encoded) == original);

        expect(!capitals_encode(std::vector<uint8_t>{'a', 0x01, 'b'}).has_value()) << "escape bytes in the text";
        expect(!capitals_encode(std::vector<uint8_t>{0x02}).has_value());
    };
    "capitals_round_trip"_test = [] {
        for (uint32_t seed = 1; seed <= 30; ++seed) {
            std::string text;
            random_text(300, 26, seed, text);
            // Uppercase some letters and add separators
            for (size_t i = 0; i < text.size(); ++i) {
                const auto r = (seed * 2654435761u + i * 40503u) >> 7;
                if (r % 5 == 0) text[i] = static_cast<char>(text[i] - 'a' + 'A');
                if (r % 7 == 0) text[i] = ' ';
                if (r % 11 == 0) text[i] = '.';
            }
            const std::vector<uint8_t> original(text.begin(), text.end());
            const auto encoded = capitals_encode(original);
            expect(encoded.has_value() && capitals_decode(*encoded) == original) << text;
        }
    };

    "bit_geometry"_test = [] {
        const std::vector<uint8_t> bytes{0b10000000, 0b00000001, 0b11000000};

        expect(stream_bit(bytes, 0) == 1 && stream_bit(bytes, 1) == 0) << "each byte from its highest bit";
        expect(stream_bit(bytes, 15) == 1 && stream_bit(bytes, 16) == 1 && stream_bit(bytes, 17) == 1);

        // At bit 17: 1 before is bit 16 (1), 2 before is bit 15 (1), 17 before is bit 0 (1), 30 before is out (0)
        expect(bit_template(bytes, 17, {1, 2, 17, 30}) == 0b1110);
        expect(bit_template(bytes, 10, {3, 10}) == 0b01);
    };

    "runs_at"_test = [] {
        // Stream bits: 00011100 11000000
        const std::vector<uint8_t> bytes{0b00011100, 0b11000000};

        // Bit 11, rows of 8: left, bit 10 (0) then bit 9 (1): a run of 1; above, bit 3 then 4 and 5 (1), 6 (0): 2
        const auto runs = runs_at(bytes, 11, 8);
        expect(runs.left == 1 && runs.above == 2);

        // Bit 6: bits 5, 4, 3 are 1, bit 2 is 0
        expect(runs_at(bytes, 6, 8).left == 3);
        expect(runs_at(bytes, 11, 8, 1).above == 1) << "capped at limit";
        expect(runs_at(bytes, 0, 8).left == 64) << "bits before the start count as 0";
    };

    "stream_bits"_test = [] {
        const std::vector<uint8_t> bytes{0b10110000, 0xFF};
        expect(stream_bits(bytes, 0, 4) == 0b1011 && stream_bits(bytes, 2, 4) == 0b1100);
        expect(stream_bits(bytes, -2, 4) == 0b0010) << "0 before the start";
        expect(stream_bits(bytes, 12, 8) == 0b11110000) << "0 beyond the end";
        expect(stream_bits(bytes, -9, 4) == 0 && stream_bits(bytes, 5, 57) >> 46 == 0b00011111111);
        expect(reverse_bits(0b1101, 4) == 0b1011 && reverse_bits(0b1101, 6) == 0b101100);
    };

    "compiled_template"_test = [] {
        std::string text;
        random_text(300, 2, 7, text);
        const std::vector<uint8_t> bytes(text.begin(), text.end());

        // Single distances, runs both ways, a run longer than one word read, distances beyond the start
        std::vector<size_t> long_run;
        for (size_t d = 1; d <= 64; ++d) {
            long_run.push_back(d);
        }
        const std::vector<std::vector<size_t>> templates{
            {1, 2, 105, 104, 103, 208}, {5}, {9, 8, 7, 6, 3, 4, 5, 300, 301}, long_run, {2000, 1999, 3000}};
        bool same = true;
        for (const auto& distances : templates) {
            const compiled_template compiled(distances);
            for (size_t pos = 0; pos <= bytes.size() * 8; pos += 7) {
                same = same && compiled.at(bytes, pos) == bit_template(bytes, pos, distances);
            }
        }
        expect(same) << "the same contexts as bit_template";
    };

    "runs_at_words"_test = [] {
        // Long runs, spanning several words, compared with a bit by bit count
        std::vector<uint8_t> bytes(64, 0);
        for (size_t i = 20; i < 40; ++i) {
            bytes[i] = 0xFF;
        }
        bytes[45] = 0b00111100;
        const auto bit = [&](const size_t pos, const size_t distance) {
            return distance <= pos ? stream_bit(bytes, pos - distance) : 0;
        };
        bool same = true;
        for (const size_t width : {2, 3, 64, 100, 200}) {
            for (const size_t limit : {1, 64, 300}) {
                for (size_t pos = 0; pos < bytes.size() * 8; pos += 3) {
                    size_t left = 1, above = 0;
                    while (left < limit && bit(pos, left) == bit(pos, left + 1)) {
                        ++left;
                    }
                    while (above < limit && above + 2 < width && bit(pos, width - above) == bit(pos, width - above - 1)) {
                        ++above;
                    }
                    const auto runs = runs_at(bytes, pos, width, limit);
                    same = same && runs.left == left && runs.above == above;
                }
            }
        }
        expect(same);
    };

    "table_bits"_test = [] {
        expect(table_bits(0) == 16) << "at least min_bits";
        expect(table_bits(3721) == 17) << "four slots per bit of data: 119,072 slots";
        expect(table_bits(131072) == 22);
        expect(table_bits(10'000'000) == 22) << "at most max_bits";
        expect(table_bits(10'000'000, 24) == 24 && table_bits(100, 22, 10) == 12);
    };

    "x86_branches"_test = [] {
        // Two calls to the target 0x1000, from offsets 0 and 0x100: different relative operands, the same bytes once absolute
        std::vector<uint8_t> code(0x110, 0x90);
        const auto call = [&](const size_t at, const uint32_t relative) {
            code[at] = 0xE8;
            for (size_t k = 0; k < 4; ++k) {
                code[at + 1 + k] = static_cast<uint8_t>(relative >> (8 * k));
            }
        };
        call(0, 0x1000);
        call(0x100, 0x1000 - 0x100);
        code[0x80] = 0xE9;                     // JMP with a far target (top byte 12): left alone
        code[0x84] = 0x12;
        const auto original = code;
        x86_branches(code, true);
        expect(std::equal(code.begin() + 1, code.begin() + 5, code.begin() + 0x101)) << "the same target, the same bytes";
        expect(code[0x81] == original[0x81] && code[0x84] == 0x12);
        x86_branches(code, false);
        expect(code == original);

        // A conditional jump (0F 85, JNE rel32) changes only on request
        code[0x40] = 0x0F;
        code[0x41] = 0x85;
        code[0x42] = 0x10;
        code[0x43] = code[0x44] = code[0x45] = 0;
        auto plain = code, with_jcc = code;
        x86_branches(plain, true);
        x86_branches(with_jcc, true, true);
        expect(plain[0x42] == 0x10 && with_jcc[0x42] == 0x50) << "0x10 + 0x40";
        x86_branches(with_jcc, false, true);
        expect(with_jcc == code);

        // Any data comes back, operands overlapping included; conditional jumps on request
        std::string text;
        random_text(3000, 26, 9, text);
        std::vector<uint8_t> noise(text.begin(), text.end());
        for (size_t i = 0; i < noise.size(); i += 7) {
            noise[i] = i % 3 == 0 ? 0xE8 : i % 3 == 1 ? 0x0F : 0x85;
            noise[(i * 31) % noise.size()] = i % 2 ? 0x00 : 0xFF;
        }
        for (const bool jcc : {false, true}) {
            auto copy = noise;
            x86_branches(copy, true, jcc);
            expect(copy != noise);
            x86_branches(copy, false, jcc);
            expect(copy == noise) << "jcc" << jcc;
        }
        std::vector<uint8_t> tiny{0xE8, 0, 0, 0};
        x86_branches(tiny, true);
        expect(tiny == std::vector<uint8_t>{0xE8, 0, 0, 0}) << "shorter than a call";

        expect(x86_branch_density(code) == 2 * 10000 / code.size()) << "two near calls, the far jump not counted";
        expect(x86_branch_density(code, 0x50, 0x20) == 0 && x86_branch_density(code, 0x200) == 0);
        expect(x86_branch_density(std::vector<uint8_t>(noise.begin(), noise.begin() + 100), 0, 0) == 0);
    };

    "x86_parser"_test = [] {
        // Hand-assembled 32-bit instructions, one of each kind the decoder tells apart
        const std::vector<std::vector<uint8_t>> instructions{
            {0x90}, {0x01, 0xC8}, {0x04, 0x05}, {0x05, 1, 2, 3, 4}, {0x06}, {0x50}, {0x74, 0x10}, {0x89, 0x45, 0xFC},
            {0x9A, 1, 2, 3, 4, 5, 6}, {0xA1, 1, 2, 3, 4}, {0xB0, 7}, {0xB8, 1, 2, 3, 4}, {0xD1, 0xE0}, {0xD9, 0xC0},
            {0xE2, 0xFE}, {0x8D, 0x04, 0x24}, {0x8B, 0x04, 0x25, 1, 2, 3, 4}, {0x68, 1, 2, 3, 4}, {0x69, 0xC0, 1, 2, 3, 4},
            {0x6A, 1}, {0x83, 0xC0, 1}, {0xC2, 8, 0}, {0xC8, 8, 0, 1}, {0xEA, 1, 2, 3, 4, 5, 6}, {0xF6, 0xC0, 1},
            {0xF7, 0xC0, 1, 2, 3, 4}, {0xF7, 0xD0}, {0xC3}, {0x66, 0xB8, 1, 2}, {0x67, 0x8B, 0x46, 1}, {0x67, 0x8B, 0x06, 1, 2},
            {0x0F, 0x05}, {0x0F, 0x84, 1, 2, 3, 4}, {0x0F, 0x70, 0xC0, 1}, {0x0F, 0xAF, 0xC1}, {0x0F, 0x38, 0x00, 0xC1},
            {0x0F, 0x3A, 0x0F, 0xC1, 8}, {0x8B, 0x80, 1, 2, 3, 4}, {0x0F, 0x1F, 0x40, 0x00}, {0xE8, 1, 2, 3, 4},
            {0x8B, 0x05, 1, 2, 3, 4}, {0xF0, 0x2E, 0x3E, 0x26, 0x36, 0x64, 0x65, 0xF2, 0xF3, 0x90}};
        const uint32_t opcodes[] = {0x90, 0x01, 0x04, 0x05, 0x06, 0x50, 0x74, 0x89, 0x9A, 0xA1, 0xB0, 0xB8, 0xD1, 0xD9, 0xE2, 0x8D,
                                    0x8B, 0x68, 0x69, 0x6A, 0x83, 0xC2, 0xC8, 0xEA, 0xF6, 0xF7, 0xF7, 0xC3, 0xB8, 0x8B, 0x8B, 0x105,
                                    0x184, 0x170, 0x1AF, 0x200, 0x30F, 0x8B, 0x11F, 0xE8, 0x8B, 0x90};
        x86_parser parser;
        bool boundaries = true;
        for (size_t i = 0; i < instructions.size(); ++i) {
            for (const auto byte : instructions[i]) {
                parser.byte_added(byte);
            }
            boundaries = boundaries && parser.next_field() == x86_parser::opcode && parser.length() == 0 && parser.previous(1) == opcodes[i];
            expect(boundaries) << "instruction" << i;
        }
        expect(parser.previous(2) == 0x8B && parser.previous(4) == 0x11F);

        // Within an instruction: MOV [EBP - 4], EAX
        parser.byte_added(0x89);
        expect(parser.next_field() == x86_parser::modrm && parser.current() == 0x89);
        parser.byte_added(0x45);
        expect(parser.next_field() == x86_parser::displacement && parser.index() == 0 && parser.modrm_byte() == 0x45);
        parser.byte_added(0xFC);
        expect(parser.previous_modrm() == 0x45);
        parser.byte_added(0x66);
        expect(parser.prefixes() == 0x66 && parser.length() == 1);

        // 64-bit: REX prefixes, 8-byte immediates and addresses
        x86_parser long_mode(true);
        const std::vector<uint8_t> code64{0x48, 0xB8, 1, 2, 3, 4, 5, 6, 7, 8, 0x48, 0x8B, 0x05, 1, 2, 3, 4,
                                          0xA1, 1, 2, 3, 4, 5, 6, 7, 8, 0x41, 0x50, 0x67, 0x8B, 0x06, 0xC3};
        for (const auto byte : code64) {
            long_mode.byte_added(byte);
        }
        expect(long_mode.next_field() == x86_parser::opcode && long_mode.previous(1) == 0xC3 && long_mode.previous(2) == 0x8B);
        expect(long_mode.previous(3) == 0x50 && long_mode.previous(4) == 0xA1);
    };

    "machine_code_detection"_test = [] {
        std::vector<uint8_t> code32, code64;
        for (int i = 0; i < 200; ++i) {
            for (const uint8_t b : {0x8B, 0x45, 0x08, 0x89, 0x04, 0x24, 0xE8, 0x10, 0x00, 0x00, 0x00}) {
                code32.push_back(b);
            }
            for (const uint8_t b : {0x48, 0x8B, 0x45, 0x08, 0x48, 0x89, 0x04, 0x24, 0x4C, 0x8D, 0x05, 0x01}) {
                code64.push_back(b);
            }
        }
        expect(!x86_long_mode(code32) && x86_long_mode(code64));
        expect(!x86_long_mode(std::vector<uint8_t>{}));

        const std::string text = "Plain text, with\ttabs and\r\nnewlines.";
        expect(control_density(std::vector<uint8_t>(text.begin(), text.end())) == 0);
        expect(control_density(std::vector<uint8_t>{0, 'a', 0x1B, 0x7F}) == 7500 && control_density({}) == 0);
        expect(control_density(std::vector<uint8_t>{'a', 'b', 0, 0}, 2) == 0) << "only the first count bytes";
    };

    "sparse_hash"_test = [] {
        const std::vector<uint8_t> history{'a', 'b', 'c', 'd'};
        expect(sparse_hash(history, 4, {2}) != sparse_hash(history, 4, {3}));
        expect(sparse_hash(history, 4, {1, 2}) != sparse_hash(history, 4, {1, 2}, 0xF0));
        expect(sparse_hash(history, 4, {2}) == sparse_hash(std::vector<uint8_t>{'x', 'y', 'c', 'z'}, 4, {2})) << "only the bytes asked for";
        expect(sparse_hash(history, 1, {2}) == sparse_hash(std::vector<uint8_t>{0, 0}, 2, {2})) << "missing bytes count as 0";
        expect(sparse_hash(history, 4, {1}, 0xF0) == sparse_hash(std::vector<uint8_t>{'a', 'b', 'c', 'o'}, 4, {1}, 0xF0)) << "masked";
    };

    "mixed_radix"_test = [] {
        // (2 * 1 + 0) * 5 + 4
        const std::vector<gmprime> digits{2, 0, 4}, radices{3, 1, 5};
        const auto packed = mixed_radix_pack(digits, radices);
        expect(packed.value == 14 && packed.radix == 15);
        expect(mixed_radix_unpack(packed.value, radices) == digits);

        const auto none = mixed_radix_pack({}, {});
        expect(none.value == 0 && none.radix == 1 && mixed_radix_unpack(0, {}).empty());

        std::vector<gmprime> many, bases;
        for (unsigned i = 0; i < 37; ++i) {
            bases.emplace_back(gmprime::binomial(40 + i, i % 9));
            many.emplace_back((bases.back() - 1) * (i % 3) / 2);
        }
        expect(mixed_radix_unpack(mixed_radix_pack(many, bases).value, bases) == many);
    };

    "rank_blocks"_test = [] {
        std::string text;
        random_text(150, 5, 11, text);
        const gmprime field = gmprime::from_bytes(std::vector<uint8_t>(text.begin(), text.end())) & ~gmprime(text.size() % 2);
        const size_t length = text.size() * 8 + 3;

        const auto blocks = rank_blocks(field, length, 64);
        expect(blocks.ones.size() == (length + 63) / 64) << "the last block shorter";
        size_t ones = 0;
        gmprime radix = 1;
        for (size_t b = 0; b < blocks.ones.size(); ++b) {
            ones += blocks.ones[b];
            radix *= gmprime::binomial(std::min<size_t>(64, length - b * 64), blocks.ones[b]);
        }
        expect(ones == field.popcount() && blocks.radix == radix && blocks.rank < radix);
        expect(unrank_blocks(blocks.ones, blocks.rank, length, 64) == field);

        expect(throws<std::invalid_argument>([&] { [[maybe_unused]] const auto r = rank_blocks(field, length, 0); }));
        expect(throws<std::invalid_argument>([&] { [[maybe_unused]] const auto r = rank_blocks(field, 8, 4); }));
    };

    "multiset_rank"_test = [] {
        const auto rank = [](const std::string& text) { return multiset_rank(std::vector<uint8_t>(text.begin(), text.end())); };
        expect(rank("aab") == 0 && rank("aba") == 1 && rank("baa") == 2);
        expect(rank("") == 0 && rank("zzz") == 0);

        std::vector<size_t> counts(256);
        counts['a'] = 2;
        counts['b'] = 1;
        expect(multiset_arrangements(counts) == 3);
        const auto text = [&](const gmprime& r) {
            const auto bytes = from_multiset_rank(r, counts);
            return std::string(bytes.begin(), bytes.end());
        };
        expect(text(0) == "aab" && text(2) == "baa");
        expect(from_multiset_rank(0, {}).empty()) << "no bytes, one arrangement";

        std::string sentence;
        random_text(200, 9, 3, sentence);
        const std::vector<uint8_t> bytes(sentence.begin(), sentence.end());
        std::vector<size_t> histogram(256);
        for (const auto byte : bytes) {
            ++histogram[byte];
        }
        const gmprime r = multiset_rank(bytes);
        expect(r < multiset_arrangements(histogram) && from_multiset_rank(r, histogram) == bytes);

        expect(throws<std::invalid_argument>([&] { [[maybe_unused]] const auto s = from_multiset_rank(3, counts); }));
        expect(throws<std::invalid_argument>([&] { [[maybe_unused]] const auto s = from_multiset_rank(-1, counts); }));
    };

    "byte_set"_test = [] {
        byte_set digits;
        for (uint8_t d = '0'; d <= '9'; ++d) {
            digits.add(d);
        }
        digits.finish();
        expect(digits.has('5') && !digits.has('a'));
        // 0x30 to 0x39: 0011 decided, then the low half open
        expect(digits.decided(1) == 0 && digits.decided(2) == 0 && digits.decided(4) == 1 && digits.decided(9) == 1);
        expect(digits.decided(19) == -1) << "the low half";
        expect(digits.decided(0b10011100) == -1) << "no digit 1100: both bits impossible";

        byte_set none;
        none.finish();
        expect(none.decided(1) == -1);
    };

    "crc32"_test = [] {
        const std::string text = "123456789";
        const std::vector<uint8_t> bytes(text.begin(), text.end());
        expect(crc32(bytes) == 0xCBF43926u) << "the standard check value";
        expect(crc32({}) == 0);

        const std::vector<uint8_t> first(bytes.begin(), bytes.begin() + 4), rest(bytes.begin() + 4, bytes.end());
        expect(crc32(rest, crc32(first)) == crc32(bytes)) << "continued over parts";
    };

    "cost_profile"_test = [] {
        cost_profile cost;
        for (int i = 0; i < 8; ++i) {
            cost.add(i % 2, 32768);
        }
        cost.add(1, 16384);

        expect(cost.per_byte().size() == 1 && std::fabs(cost.per_byte()[0] - 8) < 1e-12) << "1 bit at p = 1/2";
        expect(std::fabs(cost.total() - 10) < 1e-12) << "and the started byte: 2 bits at p = 1/4";
    };
};

namespace {

uint64_t double_bits(const double v) {
    uint64_t b;
    std::memcpy(&b, &v, 8);
    return b;
}

uint64_t float_bits(const float v) {
    uint32_t b;
    std::memcpy(&b, &v, 4);
    return b;
}

// Never throws in practice (small tests): noexcept spares its callers a cleanup path
void put_field(std::vector<uint8_t>& out, const uint64_t bits, const int size, const bool big_endian) noexcept {
    for (int i = 0; i < size; ++i) {
        out.push_back(static_cast<uint8_t>(bits >> (8 * (big_endian ? size - 1 - i : i))));
    }
}

// Records of 24 bytes: an identifier, a decimal double in little-endian (one in 50 off its lattice), one in
// big-endian, an IBM single precision value (k / 64) and a float on a lattice of 0.25
std::vector<uint8_t> lattice_records(const size_t count) {
    std::vector<uint8_t> out;
    uint32_t seed = 7;
    int64_t a = 100000, b = -3000;
    for (uint32_t r = 0; r < count; ++r) {
        seed = seed * 1664525u + 1013904223u;
        a += (seed >> 8) % 40;
        b += static_cast<int64_t>((seed >> 16) % 9) - 4;
        put_field(out, r, 4, false);
        put_field(out, double_bits(r % 50 == 3 ? 1.0 / 3 : static_cast<double>(a) * 0.001), 8, false);
        put_field(out, double_bits(static_cast<double>(b) * 0.01), 8, true);
        put_field(out, *GMPrime::coding::ibm_lattice_bits(1 + (seed >> 4) % 3000, -6), 4, true);
        put_field(out, float_bits(static_cast<float>((seed >> 12) % 200) * 0.25f), 4, false);
    }
    return out;
}

// 16-bit pixels of a smooth image, rows of width pixels
std::vector<uint8_t> image16(const size_t width, const size_t height, const bool big_endian) {
    std::vector<uint8_t> out;
    uint32_t seed = 3;
    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            seed = seed * 1664525u + 1013904223u;
            put_field(out, 1000 + 3 * x + 5 * y + (x * y) % 7 + (seed >> 20) % 4, 2, big_endian);
        }
    }
    return out;
}

} // namespace

suite<"gmprime::records"> records_suite = [] {
    using namespace GMPrime::coding;

    "record_width"_test = [] {
        std::vector<uint8_t> data;
        for (int i = 0; i < 700; ++i) {
            for (const uint8_t b : {uint8_t(i), uint8_t(i >> 8), uint8_t{7}, uint8_t{0}, uint8_t{9}, uint8_t{1}, uint8_t{2}}) {
                data.push_back(b);
            }
        }
        const auto eq = byte_autocorrelation(data, 20);
        expect(eq.size() == 21 && eq[7] > 0.7 && eq[14] > 0.7 && eq[1] < 0.3);
        expect(byte_autocorrelation(data, 20, 1 << 20, 3)[7] == eq[7] || byte_autocorrelation(data, 20, 1 << 20, 3)[7] > 0.7) << "every third byte";
        expect(byte_autocorrelation({1, 2, 3}, 10)[5] == 0) << "beyond the data";

        const auto records = find_record_width(data, 64);
        expect(records.width == 7 && records.score > records.median);
        expect(find_record_width(data, 1).width == 0) << "too few distances";
    };

    "image16"_test = [] {
        const auto little = find_image16(image16(300, 200, false));
        expect(little.found && little.little_endian && little.start == 0 && little.row_bytes == 600 && little.row_score < 0.85);
        const auto big = find_image16(image16(160, 300, true));
        expect(big.found && !big.little_endian && big.start == 0 && big.row_bytes == 320);

        // A slow ramp: samples, but no rows
        std::vector<uint8_t> ramp;
        for (size_t i = 0; i < 40000; ++i) {
            put_field(ramp, 20000 + i / 16, 2, false);
        }
        const auto line = find_image16(ramp);
        expect(line.found && line.row_bytes == 0);

        std::string text;
        random_text(40000, 26, 5, text);
        expect(!find_image16(std::vector<uint8_t>(text.begin(), text.end())).found) << "no 16-bit samples";
        expect(!find_image16(image16(10, 10, false)).found) << "too short";
    };

    "numeric_prediction"_test = [] {
        // Learns that the value is twice the feature
        lms_predictor lms(2, 3);
        expect(lms.size() == 2);
        int64_t last_error = 0;
        for (int i = 0; i < 400; ++i) {
            const int64_t x = 100 + (i * 37) % 900;
            const int64_t features[] = {x, 1};
            last_error = std::abs(lms.predict(features) - 2 * x);
            lms.update(2 * x);
        }
        expect(last_error < 40) << last_error;

        expect(quantize_difference(0) == 0 && quantize_difference(1) == 1 && quantize_difference(-1) == -1);
        expect(quantize_difference(2) == 2 && quantize_difference(3) == 3 && quantize_difference(5) == 4 && quantize_difference(6) == 5);
        expect(quantize_difference(8) == 6 && quantize_difference(-11) == -6);
        expect(quantize_difference(std::numeric_limits<int64_t>::min()) == -126 && quantize_difference(std::numeric_limits<int64_t>::max()) == 125);

        std::vector<uint8_t> pairs{1, 2, 3, 4, 5};
        swap_pairs(pairs, 0);
        expect(pairs == std::vector<uint8_t>{2, 1, 4, 3, 5});
        swap_pairs(pairs, 1);
        expect(pairs == std::vector<uint8_t>{2, 4, 1, 5, 3});

        expect(bit_middle(0b1, 1, 4) == 12 && bit_middle(0, 0, 8) == 128);
    };

    "float_fields"_test = [] {
        uint8_t bytes[] = {1, 2, 3, 4};
        expect(read_field(bytes, {4, false, false}) == 0x04030201u && read_field(bytes, {4, true, false}) == 0x01020304u);
        write_field(bytes, {4, true, false}, 0xA1B2C3D4u);
        expect(bytes[0] == 0xA1 && bytes[3] == 0xD4);
        write_field(bytes, {4, false, false}, 0xA1B2C3D4u);
        expect(bytes[0] == 0xD4 && bytes[3] == 0xA1);

        const float_format ibm{4, true, true}, single{4, false, false}, binary64{8, false, false};
        expect(field_value(0x41100000, ibm) == 1.0 && field_value(0xC1100000, ibm) == -1.0);
        expect(field_value(double_bits(2.5), binary64) == 2.5 && field_value(float_bits(0.75f), single) == 0.75);
        expect(value_bits(2.5, binary64) == double_bits(2.5) && value_bits(0.75, single) == float_bits(0.75f));
        expect(value_bits(1.0 / 64, ibm) == 0x3F400000) << "IBM steps: powers of 2";

        expect(zigzag(0) == 0 && zigzag(-1) == 1 && zigzag(1) == 2 && unzigzag(zigzag(-12345)) == -12345);
        expect(unzigzag(zigzag(std::numeric_limits<int64_t>::min())) == std::numeric_limits<int64_t>::min());
    };

    "ibm_lattice"_test = [] {
        // 1 = 1/16 * 16^1, 1/64 = 1/4 * 16^-1: normalized, first hexadecimal digit not 0
        expect(*ibm_lattice_bits(64, -6) == 0x41100000 && *ibm_lattice_bits(-64, -6) == 0xC1100000);
        expect(*ibm_lattice_bits(1, -6) == 0x3F400000 && *ibm_lattice_bits(0, -6) == 0);
        expect(!ibm_lattice_bits((1 << 25) + 1, 0)) << "more bits than the fraction holds";
        expect(!ibm_lattice_bits(1, 300) && !ibm_lattice_bits(1, -300)) << "exponent out of range";

        expect(*ibm_step_power(0x3F400000) == -6 && *ibm_step_power(0x41100000) == 0);
        expect(!ibm_step_power(0x3F300000) && !ibm_step_power(0xBF400000) && !ibm_step_power(0x3F000000));
    };

    "lattice_bits"_test = [] {
        const auto product = multiply_wide(~uint64_t{0}, ~uint64_t{0});
        expect(product.high == ~uint64_t{1} && product.low == 1);

        // k * step as the processor multiplies it, k exact in the format
        const float_format binary64{8, false, false}, single{4, false, false};
        uint32_t seed = 11;
        bool same = true;
        for (int i = 0; i < 20000; ++i) {
            seed = seed * 1664525u + 1013904223u;
            const int64_t k = static_cast<int64_t>(seed >> 4) * (i % 2 ? 1 : -1) * (i % 3 == 0 ? 1001 : 1);
            const double step = std::ldexp(1.0 + (seed % 9973) / 9973.0, static_cast<int>(seed % 50) - 30);
            same = same && *lattice_bits(k, double_bits(step), binary64) == double_bits(static_cast<double>(k) * step);
            const auto small = static_cast<int64_t>(seed >> 12);
            const float fstep = static_cast<float>(step);
            same = same && *lattice_bits(small, float_bits(fstep), single) == float_bits(static_cast<float>(small) * fstep);
        }
        expect(same);

        expect(*lattice_bits(0, double_bits(0.01), binary64) == 0);
        expect(*lattice_bits(1, double_bits(1.0), binary64) == double_bits(1.0)) << "no rounding";
        // (2^53 + 1) (2 - 2^-52) = 2^54 - 2^-52: the significand rounds up to 2^53, the exponent grows
        expect(*lattice_bits((int64_t{1} << 53) + 1, double_bits(std::nextafter(2.0, 0.0)), binary64) == double_bits(std::ldexp(1.0, 54)));
        expect(*lattice_bits(int64_t{1} << 60, 1, binary64) == double_bits(std::ldexp(1.0, -1014))) << "subnormal step";
        expect(!lattice_bits(3, 1, binary64)) << "subnormal result";
        expect(!lattice_bits(int64_t{1} << 40, double_bits(1e300), binary64)) << "overflow";
        expect(!lattice_bits(1, double_bits(std::numeric_limits<double>::infinity()), binary64) && !lattice_bits(1, 0, binary64));
        expect(*lattice_bits(64, 0x3F400000, {4, true, true}) == 0x41100000 && !lattice_bits(1, 0x3F300000, {4, true, true})) << "IBM";
    };

    "lattice_points"_test = [] {
        lattice_field field;
        field.step_bits = double_bits(0.01);
        field.residual_bits = 3;
        const auto exact = to_lattice(double_bits(1234 * 0.01), field, 0.01);
        expect(exact && exact->k == 1234 && exact->residual == 0);
        const auto near = to_lattice(double_bits(1234 * 0.01) + 2, field, 0.01);
        expect(near && near->k == 1234 && near->residual == 2);
        expect(!to_lattice(double_bits(1234 * 0.01) + 9, field, 0.01)) << "residual beyond 3 bits";
        expect(!to_lattice(double_bits(std::nan("")), field, 0.01) && !to_lattice(double_bits(1e300), field, 0.01));
        expect(!to_lattice(double_bits(-0.0), field, 0.01)) << "-0 is not k = 0";
        lattice_field broken = field;
        broken.step_bits = 0;
        expect(!to_lattice(double_bits(0.5), broken, 0.01));
        field.residual_bits = 0;
        expect(!to_lattice(double_bits(1234 * 0.01) + 1, field, 0.01) && to_lattice(double_bits(1234 * 0.01), field, 0.01));

        field.residual_bits = 4;
        for (const bool offset_binary : {false, true}) {
            field.offset_binary = offset_binary;
            for (const int64_t k : {int64_t{0}, int64_t{-5}, int64_t{123456}}) {
                const auto packed = pack_point(field, k, -3);
                expect(packed && unpack_point(field, *packed).k == k && unpack_point(field, *packed).residual == -3);
            }
            expect(!pack_point(field, int64_t{1} << 62, 0)) << "k beyond the field";
            expect(!pack_point(field, -(int64_t{1} << 62), 0));
        }
        field.residual_bits = 0;
        const int64_t extremes[] = {77, std::numeric_limits<int64_t>::min(), std::numeric_limits<int64_t>::max()};
        bool whole = true;
        for (const auto k : extremes) {
            whole = whole && unpack_point(field, *pack_point(field, k, 0)).k == k;
        }
        expect(whole) << "64 bits for k, offset binary";
    };

    "lattice_step"_test = [] {
        std::vector<double> values;
        for (int i = 1; i <= 300; ++i) {
            values.push_back(((i * 7919) % 1000) * 0.25 + 1);
        }
        const auto step = find_lattice_step(values);
        expect(step && std::fabs(*step - 0.25) < 1e-9);
        expect(!find_lattice_step(std::vector<double>(10, 1.5))) << "too few values";
        expect(!find_lattice_step(std::vector<double>(100, 1.5))) << "a single value";
        std::vector<double> scattered;
        for (int i = 1; i <= 300; ++i) {
            scattered.push_back(std::sqrt(static_cast<double>(i * 7919 % 10007)) * 1.234567);
        }
        expect(!find_lattice_step(scattered)) << "no lattice";
    };

    "decimal_steps"_test = [] {
        const float_format binary64{8, false, false};
        std::vector<uint64_t> values;
        uint32_t seed = 21;
        for (int i = 0; i < 400; ++i) {
            seed = seed * 1664525u + 1013904223u;
            // One in five off the lattice
            values.push_back(double_bits(i % 5 == 4 ? std::sqrt(static_cast<double>(seed % 1000) + 2) : static_cast<double>(seed % 100000) * 0.01));
        }
        const auto step = find_decimal_step(values, binary64);
        expect(step && *step == 0.01);
        expect(!find_decimal_step(std::vector<uint64_t>(values.begin(), values.begin() + 30), binary64)) << "too few values";
        expect(!find_decimal_step(values, {4, true, true})) << "IBM steps are powers of 2";
        std::vector<uint64_t> scattered;
        for (int i = 1; i <= 300; ++i) {
            scattered.push_back(double_bits(std::sqrt(static_cast<double>(i * 7919 % 10007)) * 1.234567));
        }
        expect(!find_decimal_step(scattered, binary64));

        std::vector<uint8_t> records(8 * 300, 7);
        for (size_t r = 0; r < 300; ++r) {
            seed = seed * 1664525u + 1013904223u;
            records[r * 8 + 3] = static_cast<uint8_t>(seed >> 24);
        }
        expect(raw_field_cost(records, 8, 0, 2, 300) < 300 && raw_field_cost(records, 8, 3, 1, 300) > 300 * 6) << "constant bytes cost little";
    };

    "lattice_window"_test = [] {
        // A header of 1000 bytes that is not records, then records of 16 bytes: two decimal doubles
        std::vector<uint8_t> data;
        uint32_t seed = 5;
        for (int i = 0; i < 1000; ++i) {
            seed = seed * 1664525u + 1013904223u;
            data.push_back(static_cast<uint8_t>(seed >> 24));
        }
        const size_t header = data.size();
        for (int r = 0; r < 1500; ++r) {
            seed = seed * 1664525u + 1013904223u;
            put_field(data, double_bits(static_cast<double>(100000 + r * 37 + seed % 30) * 1e-7), 8, false);
            put_field(data, double_bits(static_cast<double>(seed % 900000) * 1e-6), 8, false);
        }
        const std::vector<uint8_t> records(data.begin() + static_cast<std::ptrdiff_t>(header), data.end());
        const auto fields = find_lattice_fields_window(data, 16);
        expect(fields.size() == 2 && fields[0].offset == 0 && fields[1].offset == 8 && header % 16 == 8) << "the records' doubles, out of phase";
        expect(find_lattice_fields_window(records, 16).size() == 2 && find_lattice_fields_window(records, 16)[0].offset == 0);
        expect(find_lattice_fields_window(records, 3).empty() && find_lattice_fields_window(lattice_records(200), 24).empty());

        // Two floats side by side are not one double; more than 8 fields take more flag bytes
        std::vector<uint8_t> wide;
        for (int r = 0; r < 600; ++r) {
            for (int j = 0; j < 12; ++j) {
                put_field(wide, float_bits(static_cast<float>((r * (j + 3) * 31) % 500) * 0.5f), 4, false);
            }
        }
        const auto twelve = find_lattice_fields_window(wide, 48);
        bool floats = twelve.size() == 12;
        for (const auto& f : twelve) {
            floats = floats && f.format.size == 4;
        }
        expect(floats) << twelve.size();
        expect(find_lattice_fields_window(wide, 48, 20000, 10).size() == 10) << "at most max_fields";
        expect(lattice_flag_bytes(12) == 2 && lattice_flag_bytes(8) == 1 && lattice_flag_bytes(0) == 0);
        const auto encoded = lattice_encode(wide, 48, twelve);
        expect(encoded.size() == wide.size() + 600 * 2 && lattice_decode(encoded, 48, twelve) == std::optional(wide));
    };

    "scan_transform"_test = [] {
        // Decimal doubles at any place, among other bytes, the marker among them
        std::vector<uint8_t> data;
        uint32_t seed = 9;
        double x = 2.0, y = 48.0;
        for (int i = 0; i < 800; ++i) {
            seed = seed * 1664525u + 1013904223u;
            const int junk = static_cast<int>(seed % 5);
            for (int j = 0; j < junk; ++j) {
                data.push_back(static_cast<uint8_t>(seed >> (8 * j)));
            }
            x += static_cast<double>(static_cast<int>(seed % 2001) - 1000) * 1e-6;
            y += static_cast<double>(static_cast<int>((seed >> 11) % 2001) - 1000) * 1e-6;
            put_field(data, double_bits(std::round(x * 1e6) / 1e6), 8, false);
            put_field(data, double_bits(std::round(y * 1e6) / 1e6), 8, false);
        }
        const auto cfg = find_scan_config(data);
        expect(!cfg.kinds.empty() && cfg.kinds[0].format.size == 8 && field_value(cfg.kinds[0].step_bits, cfg.kinds[0].format) == 1e-6);
        bool lagged = false;
        for (const auto& k : cfg.kinds) {
            lagged = lagged || k.lag == 2;
        }
        expect(lagged) << "x and y alternate: each against the one two values before";

        for (const bool in_place : {false, true}) {
            auto c = cfg;
            c.in_place = in_place;
            const auto encoded = scan_any_encode(data, c);
            expect(scan_any_decode(encoded, c) == std::optional(data)) << "in place" << in_place;
        }
        // A residual of 40 bits leaves no room for large codes: those values stay as they are
        auto roomy = cfg;
        for (auto& k : roomy.kinds) {
            k.residual_bits = 40;
        }
        expect(scan_decode(scan_encode(data, roomy), roomy) == std::optional(data) && scan_place_decode(scan_place_encode(data, roomy), roomy) == std::optional(data));

        // Every kind of lag
        for (const int lag : {0, 1, 2}) {
            auto c = cfg;
            for (auto& k : c.kinds) {
                k.lag = lag;
            }
            expect(scan_decode(scan_encode(data, c), c) == std::optional(data) && scan_place_decode(scan_place_encode(data, c), c) == std::optional(data));
        }

        std::string text;
        random_text(4000, 26, 3, text);
        expect(find_scan_config(std::vector<uint8_t>(text.begin(), text.end())).kinds.empty()) << "no doubles in text";
        expect(find_scan_config(std::vector<uint8_t>(data.begin(), data.begin() + 500)).kinds.empty()) << "too short";

        // Invalid streams
        const uint8_t m = cfg.marker;
        expect(!scan_decode({1, m}, cfg) && !scan_decode({m, 0x0F}, cfg) && !scan_decode({m, 0x71, 1}, cfg)) << "cut, bad kind, cut code";
        expect(scan_decode({m, 0}, cfg) == std::optional(std::vector<uint8_t>{m})) << "the marker itself";
        auto broken = cfg;
        broken.kinds[0].step_bits = 0;
        expect(!scan_decode({m, 0x01, 2}, broken)) << "no value for the kind";
        const auto place = [&](std::vector<uint8_t> side) {
            std::vector<uint8_t> e(32, 0);
            e.insert(e.end(), side.begin(), side.end());
            for (int b = 3; b >= 0; --b) {
                e.push_back(static_cast<uint8_t>(side.size() >> (8 * b)));
            }
            return scan_place_decode(e, cfg);
        };
        expect(!scan_place_decode({1, 2}, cfg) && !scan_place_decode({0, 0, 0, 9}, cfg)) << "no length, side too long";
        expect(!place({0x80}) && !place({0}) && !place({0, 9, 1}) && !place({0, 1, 0x80})) << "gap cut, kind missing, bad kind, count cut";
        expect(!place({0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x01, 1, 1})) << "a gap beyond the data";
        expect(!place({30, 1, 1}) && !place({0, 1, 5})) << "values beyond the data";
        auto broken_place = cfg;
        broken_place.kinds[0].step_bits = 0;
        std::vector<uint8_t> e(32, 0);
        e[7] = 4;
        for (const uint8_t b : {uint8_t{0}, uint8_t{1}, uint8_t{1}, uint8_t{0}, uint8_t{0}, uint8_t{0}, uint8_t{3}}) {
            e.push_back(b);
        }
        expect(!scan_place_decode(e, broken_place));
    };

    "numeric_text"_test = [] {
        const auto bytes = [](const std::string& t) { return std::vector<uint8_t>(t.begin(), t.end()); };
        const std::string text = "x=-12.50, ip 192.168.0.1 v2.3.4 on 2005.06.03: 1e-5 +3 -0 007 .5 a-5 x-7 9.-1 "
                                 "1234567890123456789012 3.14159265358979323846 12345678.1234567 1.1234567 5. -";
        const auto marker = numeric_text_marker(bytes(text));
        expect(marker && *marker == 0x1A);
        const auto encoded = numeric_text_encode(bytes(text), *marker);
        expect(numeric_text_decode(encoded, *marker) == std::optional(bytes(text)));

        // -12.50: the marker, sign | (2 * 18 + 2 - 1), then the digits 1, 2, 5, 0 in their code
        const std::vector<uint8_t> first{0x1A, 0x80 | 37, digit_code(1), digit_code(2), digit_code(5), digit_code(0)};
        expect(std::equal(first.begin(), first.end(), encoded.begin() + 2)) << "after x=";

        // Random texts of digits, signs, points and letters come back
        bool all = true;
        uint32_t seed = 5;
        for (int t = 0; t < 300; ++t) {
            std::string random;
            for (int i = 0; i < 60; ++i) {
                seed = seed * 1664525u + 1013904223u;
                random += "0123456789.-+e,a \n"[(seed >> 16) % 18];
            }
            const auto e = numeric_text_encode(bytes(random), 0x1A);
            all = all && numeric_text_decode(e, 0x1A) == std::optional(bytes(random));
        }
        expect(all);

        // Invalid streams: a marker at the end, layout 126, a missing long layout, no integer digit, digits missing
        // or not digits
        const auto invalid = [](std::vector<uint8_t> e) { return !numeric_text_decode(e, 0x1A); };
        expect(invalid({'a', 0x1A}) && invalid({0x1A, 126}) && invalid({0x1A, 127, 2}) && invalid({0x1A, 127, 0, 1, digit_code(1)}));
        expect(invalid({0x1A, 1, digit_code(1)}) && invalid({0x1A, 0, 'x'}) && invalid({0x1A, 127, 10, 10}));

        std::vector<uint8_t> every(256);
        for (int b = 0; b < 256; ++b) {
            every[static_cast<size_t>(b)] = static_cast<uint8_t>(b);
        }
        expect(!numeric_text_marker(every)) << "no free byte";
        every.erase(every.begin() + 0x80);
        std::vector<uint8_t> controls(every.begin(), every.begin() + 0x80);
        expect(numeric_text_marker(controls) == std::optional<uint8_t>(0x80)) << "after the control bytes, the high ones";

        expect(digit_density(bytes("ab12")) == 5000 && digit_density({}) == 0 && digit_density(bytes("12ab"), 2) == 10000);
        bool codes = true;
        for (int d = 0; d < 10; ++d) {
            codes = codes && code_digit(digit_code(d)) == d;
        }
        expect(codes && code_digit(0x31) == -1 && code_digit('a') == -1);
        expect(split_at(1).low == 0 && split_at(1).split == 5 && split_at(1).high == 9);
        expect(split_at(2).high == 4 && split_at(2).split == 2 && split_at(3).low == 5 && split_at(3).split == 7);
        expect(power10(0) == 1 && power10(19) == 10000000000000000000ull);
    };

    "decimal_prediction"_test = [] {
        expect(rescale_decimal(1234, 2, 4) == 123400 && rescale_decimal(1235, 3, 1) == 12 && rescale_decimal(1250, 3, 1) == 13);
        expect(rescale_decimal(-1250, 3, 1) == -13) << "half away from zero";
        expect(rescale_decimal(int64_t{100000000000000000}, 0, 5) == int64_t{1000000000000000000}) << "within 10^18";
        expect(decimal_position(5, 100, 10) == 68 && decimal_position(260, 100, 10) == 83) << "below, above";
        expect(decimal_position(137, 100, 10) == 14 && decimal_position(5, 0, 1) == 20) << "digit 3, its third quarter";
    };

    "float32_print"_test = [] {
        // 127.49f prints 127.489998 with 6 decimals, 127.49 with 2
        expect(float32_print(127489998, 6) == std::optional<uint64_t>(127489998) && float32_print(12749, 2) == std::optional<uint64_t>(12749));
        expect(float32_print(127490000, 6) == std::optional<uint64_t>(127489998)) << "what the float near it prints";
        const auto one = nearest_float32(1, 18);
        expect(one.m == 9671407 && one.e == -83) << "10^-18";
        const auto tie = nearest_float32(167772155, 1);
        expect(tie.m == (1u << 23) && tie.e == 1) << "16777215.5 rounds to 2^24";
        expect(nearest_float32(0, 3).m == 0);

        // Against the C library's printing of random floats
        bool same = true;
        uint32_t seed = 9;
        char buffer[64];
        for (int i = 0; i < 20000; ++i) {
            seed = seed * 1664525u + 1013904223u;
            const uint32_t bits = (seed & 0x807FFFFF) | ((90 + (seed >> 23) % 70) << 23);
            float x;
            std::memcpy(&x, &bits, 4);
            const int f = static_cast<int>(seed % 10);
            std::snprintf(buffer, sizeof buffer, "%.*f", f, static_cast<double>(std::fabs(x)));
            std::string digits;
            for (const char* c = buffer; *c; ++c) {
                if (*c != '.') {
                    digits += *c;
                }
            }
            const uint64_t n = std::stoull(digits);
            const auto x32 = nearest_float32(n, f);
            // The print of the float nearest to the printed value is the print; with 9 digits or more, that float is x
            same = same && (n < 100000000 || std::ldexp(static_cast<double>(x32.m), x32.e) == static_cast<double>(std::fabs(x)))
                   && float32_print(n, f) == std::optional<uint64_t>(n);
        }
        expect(same);

        expect(!render_decimal({1u << 23, 0}, 13) && !render_decimal({1u << 23, 50}, 0)) << "beyond 64 bits";
        expect(!render_decimal({9094948, 40}, 0) && render_decimal({9094947, 40}, 0) == std::optional<uint64_t>(9999999980506447872ull)) << "10^19 or more";
        expect(render_decimal({1, -200}, 0) == std::optional<uint64_t>(0) && render_decimal({3, -1}, 0) == std::optional<uint64_t>(2));
        expect(!render_decimal({(1u << 24) - 1, -1}, 18) && !render_decimal({(1u << 24) - 1, -10}, 15));
        expect(render_decimal({(1u << 24) - 1, -70}, 18) == std::optional<uint64_t>(14211)) << "(2^24 - 1) 2^-70, 18 decimals";
    };

    "wide_division"_test = [] {
        bool same = true;
        uint64_t seed = 77;
        for (int i = 0; i < 20000; ++i) {
            seed = seed * 6364136223846793005ull + 1442695040888963407ull;
            const uint64_t d = (seed >> (i % 60)) | 1;
            seed = seed * 6364136223846793005ull + 1442695040888963407ull;
            const uint64_t high = seed % d, low = seed * 31 + static_cast<uint64_t>(i);
            const auto [q, r] = GMPrime::coding::detail::divide_wide(high, low, d);
            const gmprime num = (gmprime(high) << 64) + gmprime(low);
            same = same && gmprime(q) == num / gmprime(d) && gmprime(r) == num % gmprime(d);
        }
        expect(same);
    };

    "civil_time"_test = [] {
        const auto epoch = civil_from_seconds(0);
        expect(epoch.year == 1970 && epoch.month == 1 && epoch.day == 1 && epoch.hour == 0);
        const auto leap = civil_from_seconds(951782400);
        expect(leap.year == 2000 && leap.month == 2 && leap.day == 29);
        const auto before = civil_from_seconds(-1);
        expect(before.year == 1969 && before.month == 12 && before.day == 31 && before.hour == 23 && before.minute == 59 && before.second == 59);
        bool same = true;
        uint64_t seed = 3;
        for (int i = 0; i < 5000; ++i) {
            seed = seed * 6364136223846793005ull + 1442695040888963407ull;
            const auto t = static_cast<std::time_t>((seed >> 32) % 4000000000u);
            const std::tm* g = std::gmtime(&t);
            const auto c = civil_from_seconds(static_cast<int64_t>(t));
            same = same && c.year == g->tm_year + 1900 && c.month == g->tm_mon + 1 && c.day == g->tm_mday && c.hour == g->tm_hour
                   && c.minute == g->tm_min && c.second == g->tm_sec;
        }
        expect(same);
    };

    "numeric_text_reader"_test = [] {
        const std::string text = "a -12.5b 1.1234567 ";
        const auto encoded = numeric_text_encode(std::vector<uint8_t>(text.begin(), text.end()), 0x1A);
        numeric_text_reader reader(0x1A);
        expect(reader.state() == numeric_text_reader::text);
        std::vector<int> states;
        int completed = 0;
        bool checked = false;
        for (const auto b : encoded) {
            if (reader.state() == numeric_text_reader::value && reader.value_index() == 2 && !checked) {
                checked = reader.value_digits() == 3 && reader.value_prefix() == 12;
            }
            if (reader.byte_added(b)) {
                ++completed;
                if (completed == 1) {
                    const auto& n = reader.number();
                    checked = checked && n.negative && n.int_digits == 2 && n.frac_digits == 1 && n.mantissa == 125;
                }
            }
            states.push_back(reader.state());
        }
        expect(checked && completed == 2 && reader.number().frac_digits == 7 && reader.number().mantissa == 11234567);
        expect(std::count(states.begin(), states.end(), numeric_text_reader::layout_f) == 1) << "the long layout";

        // Any bytes are read: layout 126 and out of range lengths are kept within bounds
        numeric_text_reader odd(0x1A), any(0x1A);
        odd.byte_added(0x1A);
        odd.byte_added(126);
        expect(odd.state() == numeric_text_reader::value && odd.number().int_digits == 1 && odd.number().frac_digits == 7);
        for (const uint8_t b : {uint8_t{0x1A}, uint8_t{127}, uint8_t{40}, uint8_t{40}}) {
            any.byte_added(b);
        }
        expect(any.state() == numeric_text_reader::value && any.number().int_digits == 18 && any.number().frac_digits == 0);
    };

    "lattice_transform"_test = [] {
        auto data = lattice_records(600);
        data.push_back(42);                     // a byte after the last record
        const auto fields = find_lattice_fields(data, 28);
        expect(fields.size() == 4) << fields.size();
        bool formats = fields.size() == 4;
        if (formats) {
            formats = fields[0].offset == 4 && fields[0].format.size == 8 && !fields[0].format.big_endian
                      && fields[1].offset == 12 && fields[1].format.big_endian && fields[2].offset == 20 && fields[2].format.ibm
                      && fields[3].offset == 24 && fields[3].format.size == 4 && !fields[3].format.ibm;
        }
        expect(formats);

        size_t escapes = 0;
        expect(lattice_cost(data, 28, fields[0], 600, &escapes) < 600 * 64 && escapes == 12) << "1 in 50 off the lattice";

        const auto encoded = lattice_encode(data, 28, fields);
        expect(encoded.size() == data.size() + 600);
        const auto decoded = lattice_decode(encoded, 28, fields);
        expect(decoded && *decoded == data);

        auto broken = fields;
        broken[1].step_bits = 0;
        expect(!lattice_decode(encoded, 28, broken)) << "no value for these fields";
        expect(broken != fields && broken[0] == fields[0] && !(broken[1] == fields[1]));

        expect(find_lattice_fields(data, 3).empty() && find_lattice_fields(lattice_records(100), 28).empty());

        // Nine float fields: at most eight, one flag bit each
        std::vector<uint8_t> wide;
        for (int r = 0; r < 400; ++r) {
            for (int j = 0; j < 9; ++j) {
                put_field(wide, float_bits(static_cast<float>((r * (j + 3) * 31) % 500) * 0.5f), 4, false);
            }
        }
        expect(find_lattice_fields(wide, 36).size() == 8);
    };
};

suite<"gmprime::iterator"> iterator = [] {
    "range_based_for_loop_read"_test = [] {
        gmprime r(4242);

        const auto ulong_adapter = r.as_span<unsigned long>();
        expect(ulong_adapter.size() == 1);

        for (const auto &it : ulong_adapter) {
            expect(it == 4242);
        }

        const auto uchar_adapter = r.as_span<unsigned char>();
        expect(uchar_adapter.size() == 2);

        int i = 0;
        for (const auto &it : uchar_adapter) {
            switch (i) {
                case 0:
                    expect(it == 146);
                    break ;
                case 1:
                    expect(it == 16);
                    break ;
                default:
                    expect(false) << "Additional block found during iteration";
            }
            i++;
        }
    };

    "range_based_for_loop_write"_test = [] {
        gmprime r(UINT16_MAX);
        int i = 0;

        for (auto &it : r.as_span<char>()) {
            it = 'a';
            it += i;
            i++;
        }

        const auto adapter =  r.as_span<char>();
        expect(adapter.size() == 2);

        i = 0;
        for (const auto& it : adapter) {
            switch (i) {
                case 0:
                    expect(it == 'a');
                break ;
                case 1:
                    expect(it == 'b');
                break ;
                default:
                    expect(false) << "Additional block found during iteration";
            }
            i++;
        }
    };

    "range_based_for_loop_reverse"_test = [] {
        gmprime r;
        int i = 0;

        r.operator[]<char>(0) = 'a';
        r.operator[]<char>(1) = 'b';
        r.operator[]<char>(2) = 'c';

        const auto adapter = r.as_span<char>(true);
        expect(adapter.size() == 3) << "span size 3 expected, got: " << adapter.size();

        for (const auto &it : adapter) {
            switch (i) {
                case 0:
                    expect(it == 'c') << "expect 'c' at idx 0 (reverse), got: " << static_cast<char>(it);
                    break ;
                case 1:
                    expect(it == 'b') << "expect 'b' at idx 1 (reverse), got: " << static_cast<char>(it);
                    break ;
                case 2:
                    expect(it == 'a') << "expect 'a' at idx 2 (reverse), got: " << static_cast<char>(it);
                    break ;
                default:
                    expect(false) << "Additional block found during iteration";
            }
            i++;
        }
    };

    "range_based_for_loop_reverse_write"_test = [] {
        gmprime r;
        int i = 0;

        r.operator[]<int>(0) = 'a';
        r.operator[]<int>(1) = 'b';
        r.operator[]<int>(2) = 'c';

        const auto adapter = r.as_span<int>(true);
        expect(adapter.size() == 3) << "span size 3 expected, got: " << adapter.size() << " " << r.to_string(2);

        for (auto &it : adapter) {
            switch (i) {
                case 0:
                    expect(it == 'c') << "expect 'c' at idx 0 (reverse), got: " << static_cast<int>(it);
                    it = 'a';
                break ;
                case 1:
                    expect(it == 'b') << "expect 'b' at idx 1 (reverse), got: " << static_cast<int>(it);
                break ;
                case 2:
                    expect(it == 'a') << "expect 'a' at idx 2 (reverse), got: " << static_cast<int>(it);
                    it = 'c';
                break ;
                default:
                    expect(false) << "Additional block found during iteration";
            }
            i++;
        }

        i = 0;
        for (const auto &it : adapter) {
            switch (i) {
                case 0:
                    expect(it == 'a') << "expect 'a' at idx 0 (reverse), got: " << static_cast<int>(it);
                break ;
                case 1:
                    expect(it == 'b') << "expect 'b' at idx 1 (reverse), got: " << static_cast<int>(it);
                break ;
                case 2:
                    expect(it == 'c') << "expect 'c' at idx 2 (reverse), got: " << static_cast<int>(it);
                break ;
                default:
                    expect(false) << "Additional block found during iteration";
            }
            i++;
        }
    };

    "range_based_for_loop_bit"_test = [] {
        gmprime v("10", 2);
        int i = 0;
        const auto adapter = v.as_span<bool>();

        expect(adapter.size() == 2);

        for (auto bit : adapter) {
            switch (i) {
                case 0:
                    expect(bit == false);
                    break ;
                case 1:
                    expect(bit == true);
                    break ;
                default:
                    expect(false) << "Additional bit found during iteration";
            }
            i++;
        }
    };

    "range_based_for_loop_bit_write"_test = [] {
        gmprime v("101", 2);
        int i = 0;

        for (auto &bit : v.as_span<bool>()) {
            bit = !bit;
        }

        // Got 010 after reverting, which directly translates to 01 (omitting leading 0 bits)

        const auto adapter = v.as_span<bool>();
        expect(adapter.size() == 2) << "Expected size of 2 got " << adapter.size();

        for (auto bit : adapter) {
            switch (i) {
                case 0:
                    expect(bit == false);
                break ;
                case 1:
                    expect(bit == true);
                break ;
                default:
                    expect(false) << "Additional bit found during iteration";
            }
            i++;
        }
    };

    "range_based_for_loop_bit_reverse"_test = [] {
        gmprime v("10", 2);
        int i = 0;
        const auto adapter = v.as_span<bool>(true);

        expect(adapter.size() == 2);

        for (auto bit : adapter) {
            switch (i) {
                case 0:
                    expect(bit == true);
                break ;
                case 1:
                    expect(bit == false);
                break ;
                default:
                    expect(false) << "Additional bit found during iteration";
            }
            i++;
        }
    };

    "iterator_copy"_test = [] {
        gmprime r(5);
        auto it = r.as_span<bool>().begin();
        const auto copy = it;

        ++it;

        expect(*copy == true);
        expect(*it == false);
    };

    "iterator_postfix_increment"_test = [] {
        gmprime r(5);
        auto it = r.as_span<bool>().begin();
        const auto old = it++;

        expect(*old == true);
        expect(*it == false);
    };

    "iterator_postfix_decrement"_test = [] {
        gmprime r(5);
        auto it = r.as_span<bool>().begin();
        ++it;
        const auto old = it--;

        expect(*old == false);
        expect(*it == true);
    };
};

suite<"gmprime::string"> string = [] {
    "to_string"_test = [] {
        const gmprime v(6548);

        expect(v.to_string() == "6548");
        expect(v.to_string(10) == "6548");
        expect(v.to_string(16) == "1994");
        expect(v.to_string(8) == "14624");
        expect(v.to_string(2) == "1100110010100");
    };

    "ostream"_test = [] {
        const gmprime v(6548);

        std::stringstream t;
        t << gmprime::dec << v;
        expect(t.str() == "6548");
        t.str("");

        t << gmprime::hex << v;
        expect(t.str() == "1994");
        t.str("");

        t << gmprime::oct << v;
        expect(t.str() == "14624");
        t.str("");

        t << gmprime::bin << v;
        expect(t.str() == "1100110010100");
        t.str("");

        t << gmprime::dec << v[0];
        expect(t.str() == "6548");
        t.str("");

        gmprime u(v);

        t << gmprime::dec << u[0];
        expect(t.str() == "6548");
        t.str("");
    };

    "ostream_format_like_int"_test = [] {
        using flags = std::ios_base::fmtflags;
        const flags bases[] = {std::ios_base::dec, std::ios_base::oct, std::ios_base::hex};
        const flags options[] = {
            flags{}, std::ios_base::showbase, std::ios_base::uppercase, std::ios_base::showpos,
            std::ios_base::showbase | std::ios_base::uppercase,
            std::ios_base::showbase | std::ios_base::showpos | std::ios_base::uppercase,
        };
        const flags adjustments[] = {flags{}, std::ios_base::left, std::ios_base::right, std::ios_base::internal};
        const long values[] = {0, 1, 42, 255, 4096, 123456789};

        const auto check = [](const flags format, const long value) {
            std::ostringstream expected, actual;
            for (auto *s : {&expected, &actual}) {
                s->flags(format);
                s->width(14);
                s->fill('*');
            }
            expected << value;
            actual << gmprime(value);
            expect(actual.str() == expected.str()) << "gmprime" << actual.str() << "int" << expected.str();
        };

        for (const auto base : bases)
            for (const auto option : options)
                for (const auto adjustment : adjustments)
                    for (const auto value : values) {
                        check(base | option | adjustment, value);
                        // An int prints its two's complement in octal and hexadecimal: compare signs in decimal only
                        if (base == std::ios_base::dec)
                            check(base | option | adjustment, -value);
                    }
    };

    "ostream_format_beyond_int"_test = [] {
        std::ostringstream t;
        t << std::showbase << std::hex << gmprime(-255) << ' ' << std::uppercase << gmprime(-255);
        expect(t.str() == "-0xff -0XFF") << "a negative number prints as a sign and its magnitude";
        t.str("");

        t << std::nouppercase << gmprime::bin << gmprime(5) << ' ' << std::uppercase << gmprime(5) << ' ' << gmprime(0);
        expect(t.str() == "0b101 0B101 0");
        t.str("");

        t << std::noshowbase << std::nouppercase << std::hex << std::setw(10) << std::internal << std::setfill('0')
          << gmprime(-255);
        expect(t.str() == "-0000000ff");
        t.str("");

        std::wostringstream w;
        w << std::showbase << std::hex << gmprime(255);
        expect(w.str() == L"0xff") << "wide streams";
    };

    "istream"_test = [] {
        std::istringstream t("  42 -17 +5 123456789012345678901234567890");
        gmprime a, b, c, d;
        t >> a >> b >> c >> d;

        expect(!t.fail());
        expect(a == 42 && b == -17 && c == 5);
        expect(d == gmprime("123456789012345678901234567890"));
    };
    "istream_bases"_test = [] {
        gmprime a, b, c, d;

        std::istringstream hex("ff 0xFF");
        hex >> std::hex >> a >> b;
        expect(a == 255 && b == 255);

        std::istringstream oct("17");
        oct >> std::oct >> c;
        expect(c == 15);

        std::istringstream bin("101 0b11");
        bin >> gmprime::bin >> a >> b;
        expect(a == 5 && b == 3);

        // No base set: detected from the prefix, like %i
        std::istringstream detect("0x1f 017 0b101 42 -0x10 0");
        detect.unsetf(std::ios_base::basefield);
        gmprime e, f;
        detect >> a >> b >> c >> d >> e >> f;
        expect(a == 31 && b == 15 && c == 5 && d == 42 && e == -16 && f == 0);
    };
    "istream_stops_and_fails"_test = [] {
        gmprime a(1), b(2);

        std::istringstream rest("42abc");
        rest >> a;
        expect(a == 42 && !rest.fail());
        expect(rest.peek() == 'a') << "reading stops at the first character that is not a digit";

        std::istringstream invalid("abc");
        invalid >> b;
        expect(invalid.fail());
        expect(b == 2) << "value modified by a failed read";

        std::istringstream end("42");
        end >> a;
        expect(end.eof() && !end.fail());

        std::istringstream sign_only("-");
        sign_only >> b;
        expect(sign_only.fail() && b == 2);

        std::istringstream blank("   ");
        blank >> b;
        expect(blank.fail() && blank.eof() && b == 2) << "nothing but whitespace, like an int";
    };
    "stream_round_trip"_test = [] {
        const gmprime original("-123456789abcdef0123456789abcdef", 16);
        std::stringstream t;
        t << std::showbase << std::hex << original;

        gmprime read;
        t.unsetf(std::ios_base::basefield);
        t >> read;
        expect(read == original);
    };

    "ostream_bin_sticky"_test = [] {
        const gmprime a(5), b(6), c(255);
        std::stringstream t;

        t << gmprime::bin << a << ' ' << b;
        expect(t.str() == "101 110") << "binary must stay active, like std::hex";
        t.str("");

        t << 5;
        expect(t.str() == "5") << "primitives are printed in decimal";
        t.str("");

        t << std::hex << c;
        expect(t.str() == "ff") << "a standard base set afterwards cancels binary";
        t.str("");

        t << gmprime::bin << std::dec << c;
        expect(t.str() == "255");
        t.str("");

        t << gmprime::bin << gmprime::hex << c;
        expect(t.str() == "ff");
        t.str("");

        std::stringstream other;
        other << c;
        expect(other.str() == "255") << "binary is a per-stream state";
    };

    "file_read_not_a_file"_test = [] {
        // Opening a directory succeeds on Linux, then tellg() returns a meaningless size (LLONG_MAX on ext4)
        const std::string dir_path = TEST_OUTPUT_DIR "a_directory";
        std::filesystem::create_directory(dir_path);
        gmprime v(42);

        expect(throws<std::runtime_error>([&] { v.load_from_file(dir_path); }));
        expect(v == 42) << "value modified by a failed load";

        std::filesystem::remove(dir_path);
    };

    "file_write_custom_gmp_allocator"_test = [] {
        const auto out_path = TEST_OUTPUT_DIR "allocator.out.bin";
        const gmprime v(0x636261);

        mp_get_memory_functions(&gmp_default_alloc, &gmp_default_realloc, &gmp_default_free);
        mp_set_memory_functions(counting_alloc, counting_realloc, counting_free);
        gmp_live_allocations = 0;

        // Sanity check of the counting itself, so that the test cannot pass for lack of counting
        {
            gmprime grown(v);
            grown <<= 1000;
            expect(gmp_live_allocations > 0) << "GMP allocations are not counted";
        }
        const auto balance_before_save = gmp_live_allocations;

        v.save_to_file(out_path);
        mp_set_memory_functions(gmp_default_alloc, gmp_default_realloc, gmp_default_free);

        expect(balance_before_save == 0) << "GMP releases are not counted";
        expect(gmp_live_allocations == 0) << "memory allocated through GMP must be released through GMP";

        std::remove(out_path);
    };

#ifdef __linux__
    "file_write_error"_test = [] {
        // /dev/full opens fine, then every write fails with ENOSPC
        const gmprime v(0x636261);
        expect(throws<std::runtime_error>([&] { v.save_to_file("/dev/full"); }));
    };
#endif

    "file_read"_test = [] {
        gmprime v;
        expect(throws([&] {v.load_from_file("files/dummy.txt");}));

        v.load_from_file(TEST_FILES_DIR "dummy.txt");
        expect(v > 0);

        expect(v.operator[]<char>(0) == 'a');
        expect(v.operator[]<char>(1) == 'b');
        expect(v.operator[]<char>(2) == 'c');
    };

    "file_write"_test = [] {
        gmprime v;
        v.load_from_file(TEST_FILES_DIR "dummy.txt");
        expect(v > 0);

        v.save_to_file(TEST_OUTPUT_DIR "dummy.out.txt");

        struct stat buffer{};
        const bool exists = (stat (TEST_OUTPUT_DIR "dummy.out.txt", &buffer) == 0);

        expect(exists) << "Out file not written";

        gmprime t;
        t.load_from_file(TEST_OUTPUT_DIR "dummy.out.txt");

        expect(t.operator[]<char>(0) == 'a');
        expect(t.operator[]<char>(1) == 'b');
        expect(t.operator[]<char>(2) == 'c');

        std::remove(TEST_OUTPUT_DIR "dummy.out.txt");
    };
};

suite<"gmprime::math"> math = [] {
    "sqrt"_test = [] {
        const gmprime v(16);
        expect(v.sqrt() == 4);
    };
    "sqrt_negative"_test = [] {
        const gmprime v(-16);
        expect(throws<std::domain_error>([&] { static_cast<void>(v.sqrt()); }));
    };
    "pow"_test = [] {
        const gmprime v(2);
        expect(v.pow(4) == 16);
    };
    "is_divisible_by_pow2"_test = [] {
        const gmprime v(16);
        expect(v.is_divisible_by_pow2(2));
        expect(v.is_divisible_by_pow2(5) == false);
    };
    "is_divisible_by"_test = [] {
        const gmprime v(16);
        expect(v.is_divisible_by(2));
        expect(v.is_divisible_by(3) == false);
    };
    "factorial"_test = [] {
        expect(gmprime::factorial(5) == 120);
    };
    "binomial"_test = [] {
        expect(gmprime::binomial(5, 2) == 10 && gmprime::binomial(0, 0) == 1 && gmprime::binomial(3, 5) == 0);
        expect(gmprime::binomial(100, 50) == gmprime("100891344545564193334812497256"));
    };
    "combination_rank"_test = [] {
        // The i-th bit set adds C(position, i): 0b1100 is the last of the C(4, 2) = 6 arrangements
        expect(gmprime(0b1).combination_rank(4) == 0 && gmprime(0b1000).combination_rank(4) == 3);
        expect(gmprime(0b0011).combination_rank(4) == 0 && gmprime(0b1100).combination_rank(4) == 5);
        expect(gmprime(0).combination_rank(10) == 0 && gmprime(0).combination_rank(0) == 0);

        // Every field of 7 bits: the ranks of each count of ones are 0 to C(7, ones) - 1, and come back
        bool exact = true;
        std::vector<std::vector<bool>> seen(8);
        for (size_t ones = 0; ones <= 7; ++ones) {
            seen[ones].resize(gmprime::binomial(7, ones).to_primitive<size_t>());
        }
        for (unsigned v = 0; v < 128; ++v) {
            const gmprime field(v);
            const gmprime rank = field.combination_rank(7);
            const size_t ones = field.popcount();
            exact = exact && rank < gmprime::binomial(7, ones) && !seen[ones][rank.to_primitive<size_t>()];
            seen[ones][rank.to_primitive<size_t>()] = true;
            exact = exact && gmprime::from_combination_rank(rank, 7, ones) == field;
        }
        expect(exact);

        // A long sparse field, with leading zeros kept by the length
        gmprime sparse;
        for (size_t p = 3; p < 2000; p += 37 + p % 11) {
            sparse.get_value() |= mpz_class(1) << static_cast<mp_bitcnt_t>(p);
        }
        const gmprime rank = sparse.combination_rank(2500);
        expect(rank < gmprime::binomial(2500, sparse.popcount()));
        expect(gmprime::from_combination_rank(rank, 2500, sparse.popcount()) == sparse);
        expect(gmprime::from_combination_rank(0, 0, 0) == 0 && gmprime::from_combination_rank(0, 5, 0) == 0);
        expect(gmprime::from_combination_rank(0, 3, 3) == 0b111) << "the only arrangement";

        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const auto r = gmprime(0b10000).combination_rank(4); }));
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const auto f = gmprime::from_combination_rank(6, 4, 2); }));
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const auto f = gmprime::from_combination_rank(0, 4, 5); }));
        expect(throws<std::invalid_argument>([] { [[maybe_unused]] const auto f = gmprime::from_combination_rank(-1, 4, 2); }));
    };
    "next_prime_number"_test = [] {
        const gmprime v(42);
        expect(v.next_prime() == 43);
    };
    "is_perfect_power"_test = [] {
        gmprime v(16);
        expect(v.is_perfect_power());

        v = 17;
        expect(!v.is_perfect_power());
    };

    "file_round_trip_trailing_zeros"_test = [] {
        const char data[] = {'a', 'b', 'c', '\0', '\0'};
        const auto in_path = TEST_OUTPUT_DIR "trailing_zeros.in.bin";
        const auto out_path = TEST_OUTPUT_DIR "trailing_zeros.out.bin";

        {
            std::ofstream in(in_path, std::ios::binary);
            in.write(data, sizeof(data));
        }

        gmprime v;
        v.load_from_file(in_path);
        v.save_to_file(out_path);

        std::ifstream out(out_path, std::ios::binary);
        const std::string written((std::istreambuf_iterator<char>(out)), std::istreambuf_iterator<char>());

        expect(written == std::string(data, sizeof(data)));

        std::remove(in_path);
        std::remove(out_path);
    };

    "file_size_reset_on_assignment"_test = [] {
        const auto in_path = TEST_OUTPUT_DIR "trailing_zeros.in.bin";
        const auto out_path = TEST_OUTPUT_DIR "trailing_zeros.out.bin";

        {
            std::ofstream in(in_path, std::ios::binary);
            in.write("a\0\0\0", 4);
        }

        gmprime v;
        v.load_from_file(in_path);
        v = 0x62;
        v.save_to_file(out_path);

        std::ifstream out(out_path, std::ios::binary);
        const std::string written((std::istreambuf_iterator<char>(out)), std::istreambuf_iterator<char>());

        expect(written == std::string("b"));

        std::remove(in_path);
        std::remove(out_path);
    };

    "file_write_negative_saves_magnitude"_test = [] {
        const auto out_path = TEST_OUTPUT_DIR "negative.out.bin";

        const gmprime v(-0x636261);
        v.save_to_file(out_path);

        std::ifstream out(out_path, std::ios::binary);
        const std::string written((std::istreambuf_iterator<char>(out)), std::istreambuf_iterator<char>());

        expect(written == std::string("abc"));

        gmprime t;
        t.load_from_file(out_path);
        expect(t == 0x636261);

        std::remove(out_path);
    };
};

suite<"gmprime::helper"> helper = [] {
    "to_primitive"_test = [] {
        gmprime v(42);
        expect(v.to_primitive<unsigned char>() == 42);
        expect(v.to_primitive<unsigned short>() == 42);
        expect(v.to_primitive<unsigned int>() == 42);
        expect(v.to_primitive<unsigned long>() == 42);
        expect(v.to_primitive<float>() == 42.0);
        expect(v.to_primitive<double>() == 42.0);

        v = -42;

        // plain char is unsigned on some platforms (e.g. Linux on ARM)
        expect(v.to_primitive<signed char>() == -42);
        expect(v.to_primitive<short>() == -42);
        expect(v.to_primitive<int>() == -42);
        expect(v.to_primitive<long>() == -42);
        expect(v.to_primitive<float>() == -42.0);
        expect(v.to_primitive<double>() == -42.0);

        v = ULONG_LONG_MAX;

        expect(v.to_primitive<unsigned char>() == UINT8_MAX);
        expect(v.to_primitive<unsigned short>() == UINT16_MAX);
        expect(v.to_primitive<unsigned int>() == UINT32_MAX);
        // unsigned long is 32 bits on some platforms (e.g. Windows)
        expect(v.to_primitive<unsigned long long>() == UINT64_MAX);
    };

    "to_primitive_64_bit_extremes"_test = [] {
        expect(gmprime("-9223372036854775808").to_primitive<int64_t>() == std::numeric_limits<int64_t>::min());
        expect(gmprime("9223372036854775807").to_primitive<int64_t>() == std::numeric_limits<int64_t>::max());
        expect(gmprime("18446744073709551615").to_primitive<uint64_t>() == std::numeric_limits<uint64_t>::max());
    };

#ifdef __SIZEOF_INT128__
    "to_primitive_wider_than_long"_test = [] {
        const auto max = ~static_cast<unsigned __int128>(0);
        expect(gmprime("340282366920938463463374607431768211455").to_primitive<unsigned __int128>() == max);
        expect(gmprime("-170141183460469231731687303715884105728").to_primitive<__int128>() == -static_cast<__int128>(max >> 1) - 1);
        expect(gmprime(-42).to_primitive<__int128>() == -42);
    };
#endif

    "hash"_test = [] {
        const std::hash<gmprime> hash;
        const gmprime big("123456789012345678901234567890");

        expect(hash(gmprime(42)) == hash(gmprime("42")));
        expect(hash(big) == hash(gmprime("123456789012345678901234567890")));
        expect(hash(gmprime()) == hash(gmprime(0)));

        // The size of the last loaded file is not part of equality, so not of the hash either
        gmprime loaded;
        loaded.load_from_file(TEST_FILES_DIR "dummy.txt");
        const gmprime rebuilt(loaded.to_string());
        expect(loaded == rebuilt);
        expect(hash(loaded) == hash(rebuilt));

        const std::unordered_set<gmprime> set{gmprime(1), gmprime(-1), gmprime(1), big};
        expect(set.size() == 3);
        expect(set.count(gmprime("-1")) == 1);

        std::unordered_map<gmprime, int> map;
        map[gmprime(5)] = 1;
        map[gmprime("5")] += 1;
        expect(map.size() == 1);
        expect(map.at(5) == 2);
    };

    "free_functions"_test = [] {
        const gmprime a(-12), b(18), zero;

        expect(abs(a) == 12);
        expect(abs(b) == 18);
        expect(abs(zero) == 0);

        expect(gcd(a, b) == 6) << "always >= 0, like std::gcd";
        expect(gcd(zero, zero) == 0);
        expect(gcd(a, 8) == 4) << "a primitive operand converts";
        expect(lcm(a, b) == 36) << "always >= 0, like std::lcm";
        expect(lcm(a, zero) == 0);

        using std::to_string;
        expect(to_string(a) == "-12");
        expect(to_string(42) == "42") << "std::to_string stays available";
    };

    "numeric_limits"_test = [] {
        using limits = std::numeric_limits<gmprime>;

        expect(limits::is_specialized);
        expect(limits::is_integer && limits::is_signed && limits::is_exact);
        expect(!limits::is_bounded && !limits::is_modulo);
        expect(limits::radix == 2);
        expect(limits::max() == 0) << "meaningless for an unbounded type, gmprime() like the primary template";
    };

    "used_words"_test = [] {
        gmprime v(ULONG_LONG_MAX);

        expect(v.used_words() == 1);
        expect(v.used_words<long>() == 1);
        expect(v.used_words<int>() == 2);
        expect(v.used_words<short>() == 4);
        expect(v.used_words<char>() == 8);

        v += 1;

        expect(v.used_words() == 2);
        expect(v.used_words<long>() == 2);
        expect(v.used_words<int>() == 3);
        expect(v.used_words<short>() == 5);
        expect(v.used_words<char>() == 9);
    };
};

int main() {
    return 0;
}
