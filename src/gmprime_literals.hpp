#ifndef GMPRIME_LITERALS_HPP
#define GMPRIME_LITERALS_HPP

namespace GMPrime {

template<typename word_t>
struct word_index {
    uint64_t index;
    static constexpr auto word_size = sizeof(word_t);
};

template <char... Cs>
[[nodiscard]] constexpr uint64_t num() {
    static_assert(((Cs >= '0' and Cs <= '9') and ...));
    uint64_t result = 0;
    for (const char c : std::array{Cs...}) {
        result = result * static_cast<uint64_t>(10) + static_cast<uint64_t>(c - '0');
    }
    return result;
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_bw() {
    return word_index<bool>{num<Cs...>()};
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_cw() {
    return word_index<char>{num<Cs...>()};
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_ucw() {
    return word_index<unsigned char>{num<Cs...>()};
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_sw() {
    return word_index<short>{num<Cs...>()};
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_usw() {
    return word_index<unsigned short>{num<Cs...>()};
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_iw() {
    return word_index<int>{num<Cs...>()};
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_uiw() {
    return word_index<unsigned int>{num<Cs...>()};
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_lw() {
    return word_index<long>{num<Cs...>()};
}

template <char... Cs>
[[nodiscard]] constexpr auto operator""_ulw() {
    return word_index<unsigned long>{num<Cs...>()};
}

}

#endif //GMPRIME_LITERALS_HPP
