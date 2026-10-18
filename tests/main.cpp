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
        if (zero) {
            taken = true;
        }
        expect(!taken);

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
        v.save_to_file(out_path);
        mp_set_memory_functions(gmp_default_alloc, gmp_default_realloc, gmp_default_free);

        expect(gmp_live_allocations == 0) << "memory allocated through GMP must be released through GMP";

        std::remove(out_path);
    };

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
