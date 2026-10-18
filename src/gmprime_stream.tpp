#pragma once

template<class Ch, class Tr, class word_t>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::reference<word_t> &obj) {
    s << static_cast<word_t>(obj);
    return s;
}

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::ios_flags &flags) {
    if (flags.base == GMPrime::ios_base::bin) {
        // Binary stays active like any base, as long as no standard base is set: clearing basefield
        // lets a later std::oct, std::dec or std::hex cancel it (primitives print in decimal meanwhile)
        s.iword(GMPrime::ios_base::binary_index()) = 1;
        s.unsetf(std::ios_base::basefield);
    } else {
        s.iword(GMPrime::ios_base::binary_index()) = 0;
        s.setf(flags.base, std::ios_base::basefield);
    }

    return s;
}

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &operator <<(std::basic_ostream<Ch, Tr> &s, const GMPrime::gmprime &obj) {
    const auto basefield = s.flags() & std::ios_base::basefield;
    int base = 10;

    if (basefield == std::ios_base::oct)
        base = 8;
    else if (basefield == std::ios_base::hex)
        base = 16;
    else if (basefield == 0 && s.iword(GMPrime::ios_base::binary_index()) != 0)
        base = 2;

    return s << obj.to_string(base);
}
