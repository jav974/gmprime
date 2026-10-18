#pragma once

#include <cctype>
#include <istream>
#include <ostream>
#include <string>

namespace GMPrime::detail {

// Sets the base of a stream, shared by the output and input manipulators
inline void apply_base(std::ios_base &s, const ios_flags &flags) {
    if (flags.base == ios_base::bin) {
        // Binary stays active like any base, as long as no standard base is set: clearing basefield
        // lets a later std::oct, std::dec or std::hex cancel it (primitives use decimal meanwhile)
        s.iword(ios_base::binary_index()) = 1;
        s.unsetf(std::ios_base::basefield);
    } else {
        s.iword(ios_base::binary_index()) = 0;
        s.setf(flags.base, std::ios_base::basefield);
    }
}

// Base of a stream: 2, 8, 10 or 16, or 0 when no base is set (detected from the prefix on input, like %i)
inline int stream_base(std::ios_base &s) {
    const auto basefield = s.flags() & std::ios_base::basefield;

    if (basefield == std::ios_base::oct)
        return 8;
    if (basefield == std::ios_base::hex)
        return 16;
    if (basefield == std::ios_base::dec)
        return 10;

    return s.iword(ios_base::binary_index()) != 0 ? 2 : 0;
}

// Value of a digit in bases up to 36, 36 for anything else
inline int digit_value(const int c) {
    if (c >= '0' && c <= '9')
        return c - '0';
    if (c >= 'a' && c <= 'z')
        return c - 'a' + 10;
    if (c >= 'A' && c <= 'Z')
        return c - 'A' + 10;

    return 36;
}

} // GMPrime::detail

template<class Ch, class Tr, class word_t>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::reference<word_t> &obj) {
    s << static_cast<word_t>(obj);
    return s;
}

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::ios_flags &flags) {
    GMPrime::detail::apply_base(s, flags);
    return s;
}

template<class Ch, class Tr>
std::basic_istream<Ch, Tr> &operator >>(std::basic_istream<Ch, Tr> &s, const GMPrime::ios_flags &flags) {
    GMPrime::detail::apply_base(s, flags);
    return s;
}

// Formats like an int: showbase (0x, 0, 0b), uppercase, showpos, and width with left, right or internal
// adjustment. A negative number prints as a sign and its magnitude in every base (-0xff), where an int
// prints its two's complement in hexadecimal and octal.
template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::gmprime &obj) {
    const auto flags = s.flags();
    const auto stream_base = GMPrime::detail::stream_base(s);
    const auto base = stream_base == 0 ? 10 : stream_base;
    const bool uppercase = (flags & std::ios_base::uppercase) != 0;

    std::string digits = obj.to_string(base);
    std::string head;

    if (digits.front() == '-') {
        head = "-";
        digits.erase(0, 1);
    } else if ((flags & std::ios_base::showpos) && base == 10) {
        // Like an int, whose octal and hexadecimal forms are unsigned
        head = "+";
    }

    // Like printf's # flag: no prefix for zero. Internal padding goes after 0x (and 0b), but before the
    // octal 0, which counts as a digit.
    if ((flags & std::ios_base::showbase) && digits != "0") {
        if (base == 16)
            head += uppercase ? "0X" : "0x";
        else if (base == 8)
            digits.insert(0, "0");
        else if (base == 2)
            head += uppercase ? "0B" : "0b";
    }

    if (uppercase) {
        for (auto &c : digits) {
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
    }

    std::basic_string<Ch, Tr> text;
    const auto length = head.size() + digits.size();
    const auto width = static_cast<size_t>(std::max<std::streamsize>(s.width(), 0));
    const std::basic_string<Ch, Tr> padding(width > length ? width - length : 0, s.fill());
    const auto adjust = flags & std::ios_base::adjustfield;

    if (adjust != std::ios_base::left && adjust != std::ios_base::internal)
        text += padding;
    for (const char c : head)
        text += s.widen(c);
    if (adjust == std::ios_base::internal)
        text += padding;
    for (const char c : digits)
        text += s.widen(c);
    if (adjust == std::ios_base::left)
        text += padding;

    s.width(0);
    return s << text;
}

// Reads like an int: optional sign, then digits in the stream's base. With no base set (unsetf(basefield)),
// the base comes from the prefix: 0x hexadecimal, 0b binary, 0 octal, else decimal. A 0x prefix is also
// accepted in hexadecimal and 0b in binary. Sets failbit and leaves the value unchanged when no digit is read.
template<class Ch, class Tr>
std::basic_istream<Ch, Tr> &operator >>(std::basic_istream<Ch, Tr> &s, GMPrime::gmprime &obj) {
    const typename std::basic_istream<Ch, Tr>::sentry sentry(s);
    if (!sentry) {
        return s;
    }

    auto *const buffer = s.rdbuf();
    int base = GMPrime::detail::stream_base(s);
    std::ios_base::iostate state = std::ios_base::goodbit;
    std::string text;
    bool has_digits = false;

    // Next character narrowed to char as an unsigned value (0xFF must not look like the end), or end
    constexpr int end = -1;
    const auto peek = [&]() -> int {
        const auto c = buffer->sgetc();
        return Tr::eq_int_type(c, Tr::eof()) ? end : static_cast<unsigned char>(s.narrow(Tr::to_char_type(c), '\0'));
    };

    int c = peek();
    if (c == '+' || c == '-') {
        if (c == '-') {
            text += '-';
        }
        buffer->sbumpc();
        c = peek();
    }

    if (c == '0' && (base == 0 || base == 16 || base == 2)) {
        buffer->sbumpc();
        c = peek();

        if ((c == 'x' || c == 'X') && (base == 0 || base == 16)) {
            base = 16;
            buffer->sbumpc();
            c = peek();
        } else if ((c == 'b' || c == 'B') && (base == 0 || base == 2)) {
            base = 2;
            buffer->sbumpc();
            c = peek();
        } else {
            // The 0 was a digit, the base octal if it was to be detected
            text += '0';
            has_digits = true;
            if (base == 0) {
                base = 8;
            }
        }
    }

    if (base == 0) {
        base = 10;
    }

    while (c != end && GMPrime::detail::digit_value(c) < base) {
        text += static_cast<char>(c);
        has_digits = true;
        buffer->sbumpc();
        c = peek();
    }

    if (c == end) {
        state |= std::ios_base::eofbit;
    }

    if (has_digits) {
        obj = GMPrime::gmprime(text, base);
    } else {
        state |= std::ios_base::failbit;
    }

    s.setstate(state);
    return s;
}

// Printed as a gmprime, so that any width works and gmprime::bin applies
template<class Ch, class Tr, bool is_const>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::bit_window<is_const> &window) {
    return s << static_cast<GMPrime::gmprime>(window);
}
