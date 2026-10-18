#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <boost/ut.hpp>
#include <gmprime.hpp>

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
        expect(gmprime(0b1101).histogram(2, 1) == std::vector<size_t>{0, 1, 1, 1});
        expect(gmprime(0b1101).histogram(2) == std::vector<size_t>{0, 1, 0, 1}) << "aligned: bits 0-1, then 2-3";
        expect(gmprime(0x0F0F).histogram(4, 8) == std::vector<size_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2})
            << "step wider than width: the low nibble of each byte";
        expect(gmprime(1).histogram(8, 0, 24)[0] == 2) << "explicit length: the zero bytes above count";
        expect(gmprime(-0b1101).histogram(2) == gmprime(0b1101).histogram(2)) << "on the magnitude";
        expect(gmprime(0b1101).histogram(8) == std::vector<size_t>(256, 0)) << "no complete window";
    };
    "histogram_invalid_width"_test = [] {
        expect(throws<std::invalid_argument>([] { static_cast<void>(gmprime(1).histogram(0)); }));
        constexpr auto too_wide = gmprime::max_histogram_width + 1;
        expect(throws<std::invalid_argument>([] { static_cast<void>(gmprime(1).histogram(too_wide)); }));
    };
    "histogram_matches_windows"_test = [&] {
        for (const auto& value : irregular) {
            for (const size_t width : {1, 3, 8, 13}) {
                for (const size_t step : {size_t{1}, width, size_t{5}}) {
                    std::vector<size_t> expected(size_t{1} << width, 0);
                    for (size_t pos = 0; pos + width <= value.used_bits(); pos += step) {
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
