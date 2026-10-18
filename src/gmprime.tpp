#pragma once

static_assert(GMP_NAIL_BITS == 0, "GMP builds with nail bits are not supported");

template<typename word_t, std::enable_if_t<std::is_arithmetic_v<word_t> && !std::is_same_v<word_t, bool>, bool>>
gmprime::gmprime(const word_t value) {
    if constexpr (std::is_floating_point_v<word_t>) {
        // GMP raises SIGFPE on infinity and NaN
        if (!std::isfinite(value)) {
            throw std::domain_error("Cannot construct from an infinite or NaN value");
        }
        mpz_set_d(_value.get_mpz_t(), static_cast<double>(value));
    } else if constexpr (sizeof(word_t) <= sizeof(long int)) {
        if constexpr (std::is_signed_v<word_t>) {
            mpz_set_si(_value.get_mpz_t(), static_cast<signed long int>(value));
        } else {
            mpz_set_ui(_value.get_mpz_t(), static_cast<unsigned long int>(value));
        }
    } else {
        // Wider than long (e.g. 64-bit integers where long is 32 bits, like Windows): write the magnitude
        using unsigned_t = std::make_unsigned_t<word_t>;
        bool negative = false;

        if constexpr (std::is_signed_v<word_t>) {
            negative = value < 0;
        }

        const auto magnitude = negative ? static_cast<unsigned_t>(0 - static_cast<unsigned_t>(value)) : static_cast<unsigned_t>(value);
        write_bits<unsigned_t>(0, sizeof(unsigned_t) * CHAR_BIT, magnitude);

        if (negative) {
            mpz_neg(_value.get_mpz_t(), _value.get_mpz_t());
        }
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

template<typename word_t>
reference<word_t> gmprime::operator[](const word_index<word_t>& index) {
    return operator[]<word_t, true>(index.index);
}

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t>, bool>>
word_adapter<word_t> gmprime::as_span(const bool reverse) {
    return word_adapter<word_t>(*this, reverse);
}

template<typename shift_t>
uintmax_t gmprime::shift_count(const shift_t shift) {
    if constexpr (std::is_signed_v<shift_t>) {
        if (shift < 0) {
            throw std::domain_error("Negative shift count");
        }
    }

    // Wider than uintmax_t (e.g. __int128): saturate, the count is out of reach anyway
    if constexpr (sizeof(shift_t) > sizeof(uintmax_t)) {
        if (shift > static_cast<shift_t>(std::numeric_limits<uintmax_t>::max())) {
            return std::numeric_limits<uintmax_t>::max();
        }
    }

    return static_cast<uintmax_t>(shift);
}

template<typename shift_t, std::enable_if_t<std::is_integral_v<shift_t> && !std::is_same_v<shift_t, bool>, bool>>
gmprime& gmprime::operator<<=(const shift_t shift) {
    shift_left(shift_count(shift));
    return *this;
}

template<typename shift_t, std::enable_if_t<std::is_integral_v<shift_t> && !std::is_same_v<shift_t, bool>, bool>>
gmprime& gmprime::operator>>=(const shift_t shift) {
    shift_right(shift_count(shift));
    return *this;
}

template<typename shift_t, std::enable_if_t<std::is_integral_v<shift_t> && !std::is_same_v<shift_t, bool>, bool>>
gmprime gmprime::operator<<(const shift_t shift) const {
    return gmprime(*this) <<= shift;
}

template<typename shift_t, std::enable_if_t<std::is_integral_v<shift_t> && !std::is_same_v<shift_t, bool>, bool>>
gmprime gmprime::operator>>(const shift_t shift) const {
    return gmprime(*this) >>= shift;
}

template<typename word_t>
size_t gmprime::used_words() const {
    constexpr auto word_size = sizeof(word_t) * CHAR_BIT;
    const auto size = used_bits();

    return size / word_size + (size % word_size != 0 ? 1 : 0);
}

template<typename word_t, std::enable_if_t<std::is_arithmetic_v<word_t> && !std::is_same_v<word_t, bool>, bool>>
word_t gmprime::to_primitive() const {
    if constexpr (std::is_floating_point_v<word_t>) {
        return static_cast<word_t>(_value.get_d());
    } else if constexpr (sizeof(word_t) <= sizeof(long int)) {
        if constexpr (std::is_signed_v<word_t>) {
            return static_cast<word_t>(_value.get_si());
        } else {
            return static_cast<word_t>(_value.get_ui());
        }
    } else {
        // Wider than long: same semantics as get_si/get_ui, low bits of the magnitude, with the sign if signed
        using unsigned_t = std::make_unsigned_t<word_t>;
        const auto magnitude = read_bits<unsigned_t>(0, sizeof(unsigned_t) * CHAR_BIT);

        if constexpr (std::is_signed_v<word_t>) {
            if (mpz_sgn(_value.get_mpz_t()) < 0) {
                return static_cast<word_t>(0 - magnitude);
            }
        }

        return static_cast<word_t>(magnitude);
    }
}

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool>>
void gmprime::set_partial_limb(const uint64_t index, const word_t value) {
    using unsigned_t = std::make_unsigned_t<word_t>;
    constexpr auto word_bits = sizeof(word_t) * CHAR_BIT;

    write_bits<unsigned_t>(index * word_bits, word_bits, static_cast<unsigned_t>(value));
}

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t> && !std::is_same_v<word_t, bool>, bool>>
word_t gmprime::get_partial_limb(const uint64_t index) const {
    using unsigned_t = std::make_unsigned_t<word_t>;
    constexpr auto word_bits = sizeof(word_t) * CHAR_BIT;

    return static_cast<word_t>(read_bits<unsigned_t>(index * word_bits, word_bits));
}

template<typename unsigned_t>
unsigned_t gmprime::read_bits(const size_t pos, const size_t width) const {
    constexpr size_t limb_bits = GMP_NUMB_BITS;
    unsigned_t result = 0;

    for (size_t done = 0; done < width;) {
        const auto offset = (pos + done) % limb_bits;
        const auto count = std::min(width - done, limb_bits - offset);
        // mpz_getlimbn reads the magnitude and returns 0 beyond the used limbs
        auto bits = mpz_getlimbn(_value.get_mpz_t(), static_cast<mp_size_t>((pos + done) / limb_bits)) >> offset;

        if (count < limb_bits) {
            bits &= (static_cast<mp_limb_t>(1) << count) - 1;
        }

        result |= static_cast<unsigned_t>(static_cast<unsigned_t>(bits) << done);
        done += count;
    }

    return result;
}

template<typename unsigned_t>
void gmprime::write_bits(const size_t pos, const size_t width, const unsigned_t value) {
    constexpr size_t limb_bits = GMP_NUMB_BITS;

    if (width == 0) {
        return;
    }

    const auto sign = mpz_sgn(_value.get_mpz_t()) < 0 ? -1 : 1;
    auto current_size = mpz_size(_value.get_mpz_t());
    const auto new_size = std::max(current_size, (pos + width - 1) / limb_bits + 1);
    const auto limbs = mpz_limbs_modify(_value.get_mpz_t(), static_cast<mp_size_t>(new_size));

    // mpz_limbs_modify does not initialize the limbs it adds
    while (current_size < new_size) {
        limbs[current_size++] = 0;
    }

    for (size_t done = 0; done < width;) {
        const auto index = (pos + done) / limb_bits;
        const auto offset = (pos + done) % limb_bits;
        const auto count = std::min(width - done, limb_bits - offset);
        const auto mask = (count < limb_bits ? (static_cast<mp_limb_t>(1) << count) - 1 : ~static_cast<mp_limb_t>(0)) << offset;
        const auto bits = static_cast<mp_limb_t>(value >> done) << offset;

        limbs[index] = (limbs[index] & ~mask) | (bits & mask);
        done += count;
    }

    // Normalizes away the high zero limbs
    mpz_limbs_finish(_value.get_mpz_t(), static_cast<mp_size_t>(new_size) * sign);
}
