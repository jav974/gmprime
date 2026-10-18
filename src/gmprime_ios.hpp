#ifndef GMPRIME_IOS_HPP
#define GMPRIME_IOS_HPP

#include <ios>

namespace GMPrime {
    struct ios_flags {
        std::ios_base::fmtflags base;
    };

    struct ios_base {
        // Never set on a stream, only identifies the binary manipulator: the binary state is kept in the
        // stream's own storage (see binary_index) instead of borrowing an implementation-defined fmtflags bit
        static constexpr auto bin = std::ios_base::basefield;
        static constexpr auto oct = std::ios_base::oct;
        static constexpr auto dec = std::ios_base::dec;
        static constexpr auto hex = std::ios_base::hex;

        // Index of the per-stream binary state in iword()
        static int binary_index() {
            static const int index = std::ios_base::xalloc();
            return index;
        }
    };
}

#endif //GMPRIME_IOS_HPP
