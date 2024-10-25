#ifndef GMPRIME_IOS_HPP
#define GMPRIME_IOS_HPP

#include <ios>

namespace GMPrime {
    struct ios_flags {
        std::ios_base::fmtflags base;
    };

    struct ios_base {
        static constexpr auto bin = static_cast<std::ios_base::fmtflags>(1L << 24);
        static constexpr auto oct = std::ios_base::oct;
        static constexpr auto dec = std::ios_base::dec;
        static constexpr auto hex = std::ios_base::hex;
    };
}

#endif //GMPRIME_IOS_HPP
