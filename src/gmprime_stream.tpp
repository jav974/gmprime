#pragma once

template<class Ch, class Tr, class word_t>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::reference<word_t> &obj) {
    s << static_cast<word_t>(obj);
    return s;
}

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::ios_flags &flags) {
    if (flags.base == GMPrime::ios_base::bin) {
        s.setf(flags.base);
    } else {
        s.setf(flags.base, std::ios_base::basefield);
    }

    return s;
}

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::gmprime &obj) {
    const auto flags = s.flags();
    int base = 10;

    if (flags & GMPrime::ios_base::bin)
        base = 2;
    else if (flags & GMPrime::ios_base::oct)
        base = 8;
    else if (flags & GMPrime::ios_base::dec)
        base = 10;
    else if (flags & GMPrime::ios_base::hex)
        base = 16;

    s << obj.to_string(base);

    // Restore previous mask
    if (base == 2) {
        s.unsetf(GMPrime::ios_base::bin);
    }

    return s;
}
