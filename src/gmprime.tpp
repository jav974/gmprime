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

namespace detail {

// Compares a with the floating point value x exactly, x not being NaN: returns < 0, 0 or > 0.
// Converting x to gmprime would truncate it, converting a to floating point would round it, so x is
// decomposed exactly as mantissa * 2^exponent with an integer mantissa (works for long double too,
// unlike mpz_cmp_d).
template<typename F>
int compare_floating(const gmprime& a, const F x) {
    if (std::isinf(x)) {
        return x > 0 ? -1 : 1;
    }

    int exponent = 0;
    // |x| = fraction * 2^exponent, with fraction in [0.5, 1) (or 0)
    F fraction = std::frexp(std::fabs(x), &exponent);
    gmprime mantissa;

    // Moves the bits of the fraction into the mantissa 32 at a time: each step is exact,
    // and the fraction runs out after at most digits / 32 + 1 steps
    while (fraction != 0) {
        fraction = std::ldexp(fraction, 32);
        const auto chunk = static_cast<uint32_t>(fraction);
        fraction -= static_cast<F>(chunk);
        mantissa = (mantissa << 32) + chunk;
        exponent -= 32;
    }

    if (std::signbit(x)) {
        mantissa = -mantissa;
    }

    if (exponent >= 0) {
        return mpz_cmp(a.get_value().get_mpz_t(), (mantissa << exponent).get_value().get_mpz_t());
    }

    return mpz_cmp((a << -exponent).get_value().get_mpz_t(), mantissa.get_value().get_mpz_t());
}

} // detail

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator==(const gmprime& lhs, const F rhs) {
    return !std::isnan(rhs) && detail::compare_floating(lhs, rhs) == 0;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator!=(const gmprime& lhs, const F rhs) {
    return !(lhs == rhs);
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator<(const gmprime& lhs, const F rhs) {
    return !std::isnan(rhs) && detail::compare_floating(lhs, rhs) < 0;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator<=(const gmprime& lhs, const F rhs) {
    return !std::isnan(rhs) && detail::compare_floating(lhs, rhs) <= 0;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator>(const gmprime& lhs, const F rhs) {
    return !std::isnan(rhs) && detail::compare_floating(lhs, rhs) > 0;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator>=(const gmprime& lhs, const F rhs) {
    return !std::isnan(rhs) && detail::compare_floating(lhs, rhs) >= 0;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator==(const F lhs, const gmprime& rhs) {
    return rhs == lhs;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator!=(const F lhs, const gmprime& rhs) {
    return rhs != lhs;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator<(const F lhs, const gmprime& rhs) {
    return rhs > lhs;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator<=(const F lhs, const gmprime& rhs) {
    return rhs >= lhs;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator>(const F lhs, const gmprime& rhs) {
    return rhs < lhs;
}

template<typename F, std::enable_if_t<std::is_floating_point_v<F>, bool> = true>
bool operator>=(const F lhs, const gmprime& rhs) {
    return rhs <= lhs;
}
