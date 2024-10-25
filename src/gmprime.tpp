#pragma once

inline gmprime::gmprime(): _value(0) {
}

template<typename word_t, std::enable_if_t<std::is_arithmetic_v<word_t> && !std::is_same_v<word_t, bool>, bool>>
gmprime::gmprime(word_t value) {
    if (std::is_floating_point_v<word_t>) {
        mpz_set_d(_value.get_mpz_t(), value);
    } else if (std::is_signed_v<word_t>) {
        mpz_set_si(_value.get_mpz_t(), static_cast<signed long int>(value));
    } else {
        mpz_set_ui(_value.get_mpz_t(), static_cast<unsigned long int>(value));
    }
}

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t>, bool>>
word_t gmprime::operator[](const long int index) const {
    return get_partial_limb<word_t>(index);
}

template<>
inline bool gmprime::operator[]<bool, true>(const long int index) const {
    return get_bit(index);
}

template<typename word_t>
word_t gmprime::operator[](const word_index<word_t>& index) const {
    return operator[]<word_t, true>(index.index);
}

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t>, bool>>
reference<word_t> gmprime::operator[](const long int index) {
    return reference<word_t>(this, index);
}

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t>, bool>>
word_adapter<word_t> gmprime::as_span(const bool reverse) {
    return word_adapter<word_t>(*this, reverse);
}

template<typename word_t>
size_t gmprime::used_words() const {
    constexpr auto word_size = sizeof(word_t) * CHAR_BIT;
    const auto size = used_bits();

    return size / word_size + (size % word_size != 0 ? 1 : 0);
}

template<typename word_t, std::enable_if_t<std::is_arithmetic_v<word_t> && !std::is_same_v<word_t, bool>, bool>>
word_t gmprime::to_primitive() const {
    if (std::is_floating_point_v<word_t>) {
        return _value.get_d();
    }

    if (std::is_signed_v<word_t>) {
        return _value.get_si();
    }

    return _value.get_ui();
}

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool>>
void gmprime::set_partial_limb(const uint64_t index, const word_t value) {
    constexpr auto word_size = sizeof(word_t);
    constexpr auto limb_size = sizeof(mp_limb_t);
    auto current_size = mpz_size(_value.get_mpz_t());
    const auto limbs_required = (index + 1) * word_size / limb_size + ((index + 1) * word_size % limb_size > 0 ? 1 : 0);
    const auto new_size = std::max(std::max(current_size, static_cast<unsigned long>(1)), limbs_required);
    const auto limbs_ptr = mpz_limbs_modify(_value.get_mpz_t(), new_size);
    const auto ptr = reinterpret_cast<word_t*>(limbs_ptr);
    const auto sign = mpz_sgn(_value.get_mpz_t()) < 0 ? -1 : 1;

    // With -fsanitize=address compiler option, new limbs contain random bits, so let's 0'em out
    while (current_size < new_size) {
        limbs_ptr[current_size++] = 0;
    }

    ptr[index] = value;

    mpz_limbs_finish(_value.get_mpz_t(), new_size * sign);
}

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool>>
word_t gmprime::get_partial_limb(const uint64_t index) const {
    constexpr auto word_size = sizeof(word_t);
    constexpr auto limb_size = sizeof(mp_limb_t);
    const auto limb = mpz_getlimbn(_value.get_mpz_t(), index * word_size / limb_size);

    if (word_size == limb_size) {
        return static_cast<word_t>(limb);
    }

    constexpr word_t mask = ~static_cast<word_t>(0);
    const auto limb_sector = index * word_size % limb_size;

    return (limb >> (limb_sector * CHAR_BIT)) & mask;
}
