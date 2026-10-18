// ReSharper disable CppNonExplicitConversionOperator
#pragma once

namespace detail {

// Integers a bit window converts to and from: unsigned ones, except bool and character types
// (not numbers, and the signedness of char and wchar_t depends on the platform)
template<typename U>
constexpr bool is_window_word_v = std::is_integral_v<U> && std::is_unsigned_v<U>
    && !std::is_same_v<U, bool> && !std::is_same_v<U, char> && !std::is_same_v<U, wchar_t>
    && !std::is_same_v<U, char16_t> && !std::is_same_v<U, char32_t>;

} // detail

template<bool is_const>
class bit_window {
    using num_t = std::conditional_t<is_const, const gmprime, gmprime>;

public:
    bit_window(num_t& num, const size_t pos, const size_t width): _num(&num), _pos(pos), _width(width) {
    }

    bit_window(const bit_window& other) = default;

    [[nodiscard]] size_t position() const {
        return _pos;
    }

    [[nodiscard]] size_t width() const {
        return _width;
    }

    // Throws std::overflow_error if the window is wider than U, whatever the value it holds
    template<typename U, std::enable_if_t<detail::is_window_word_v<U>, bool> = true>
    operator U() const {
        constexpr auto digits = static_cast<size_t>(std::numeric_limits<U>::digits);

        if (_width > digits) {
            throw std::overflow_error("A " + std::to_string(_width) + "-bit window does not fit in a "
                                      + std::to_string(digits) + "-bit integer");
        }

        return _num->template read_bits<U>(_pos, _width);
    }

    operator gmprime() const {
        return _num->read_window(_pos, _width);
    }

    // Keeps the low width bits of value, zero-extended if the window is wider
    template<typename U, bool C = is_const, std::enable_if_t<detail::is_window_word_v<U> && !C, bool> = true>
    bit_window& operator=(const U value) {
        if (_width <= static_cast<size_t>(std::numeric_limits<U>::digits)) {
            _num->template write_bits<U>(_pos, _width, value);
        } else {
            _num->write_window(_pos, _width, gmprime(value));
        }

        return *this;
    }

    // Keeps the low width bits of the magnitude of value
    template<bool C = is_const, std::enable_if_t<!C, bool> = true>
    bit_window& operator=(const gmprime& value) {
        _num->write_window(_pos, _width, value);
        return *this;
    }

    // Copies the bits of the other window, like std::bitset::reference: assigning never rebinds a window
    bit_window& operator=(const bit_window& other) {
        static_assert(!is_const, "a window on a const gmprime is read-only");
        return *this = static_cast<gmprime>(other);
    }

    // Signed integers, floating point, bool and character types are refused: a sign has no meaning
    // within a window (without this, they would silently convert to gmprime and lose their sign)
    template<typename S, std::enable_if_t<std::is_arithmetic_v<S> && !detail::is_window_word_v<S>, bool> = true>
    bit_window& operator=(S value) = delete;

private:
    num_t* _num;
    size_t _pos;
    size_t _width;
};
