#include <iostream>
#include <sys/stat.h>
#include <boost/ut.hpp>
#include <gmprime.hpp>

using namespace boost::ut;
using namespace GMPrime;

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
    "unary+"_test = [] {
        const gmprime r(-42);
        expect(+r == -42);
    };
    "unary-"_test = [] {
        const gmprime r(42);
        expect(-r == -42);
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

    "file_read"_test = [] {
        gmprime v;
        expect(throws([&] {v.load_from_file("files/dummy.txt");}));

        v.load_from_file("../tests/files/dummy.txt");
        expect(v > 0);

        expect(v.operator[]<char>(0) == 'a');
        expect(v.operator[]<char>(1) == 'b');
        expect(v.operator[]<char>(2) == 'c');
    };

    "file_write"_test = [] {
        gmprime v;
        v.load_from_file("../tests/files/dummy.txt");
        expect(v > 0);

        v.save_to_file("../tests/files/dummy.out.txt");

        struct stat buffer{};
        const bool exists = (stat ("../tests/files/dummy.out.txt", &buffer) == 0);

        expect(exists) << "Out file not written";

        gmprime t;
        t.load_from_file("../tests/files/dummy.out.txt");

        expect(t.operator[]<char>(0) == 'a');
        expect(t.operator[]<char>(1) == 'b');
        expect(t.operator[]<char>(2) == 'c');

        std::remove("../tests/files/dummy.out.txt");
    };
};

suite<"gmprime::math"> math = [] {
    "sqrt"_test = [] {
        const gmprime v(16);
        expect(v.sqrt() == 4);
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

        expect(v.to_primitive<char>() == -42);
        expect(v.to_primitive<short>() == -42);
        expect(v.to_primitive<int>() == -42);
        expect(v.to_primitive<long>() == -42);
        expect(v.to_primitive<float>() == -42.0);
        expect(v.to_primitive<double>() == -42.0);

        v = ULONG_LONG_MAX;

        expect(v.to_primitive<unsigned char>() == UINT8_MAX);
        expect(v.to_primitive<unsigned short>() == UINT16_MAX);
        expect(v.to_primitive<unsigned int>() == UINT32_MAX);
        expect(v.to_primitive<unsigned long>() == UINT64_MAX);
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
