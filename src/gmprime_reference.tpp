// ReSharper disable CppNonExplicitConversionOperator
#pragma once

template<typename T, std::enable_if_t<std::is_integral_v<T>, bool>>
class reference {
public:
    reference() = default;
    reference(const reference &other) = default;
    reference(reference &&other) noexcept = default;
    explicit reference(gmprime* const num, const uint64_t index): _num(num), _index(index) {
    }

    template<typename U = T, std::enable_if_t<!std::is_same_v<U, bool>, bool> = true>
    operator U() const {
        return _num->get_partial_limb<U>(_index);
    }

    operator bool() const {
        return _num->get_bit(_index);
    }

    template<typename U = T, std::enable_if_t<std::is_integral_v<U> && !std::is_same_v<U, bool>, bool> = true>
    reference& operator=(const U value) {
        _num->set_partial_limb<T>(_index, static_cast<T>(value));
        return *this;
    }

    reference& operator=(const bool value) {
        _num->set_bit(_index, value);
        return *this;
    }

    reference& operator=(const reference& other) {
        if (this != &other) {
            _num = other._num;
            _index = other._index;
        }
        return *this;
    }

    template<typename U = T, std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
    bool operator==(const U value) const {
        return operator T() == value;
    }

    template<typename U = T, std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
    bool operator!=(const U value) const {
        return operator T() != value;
    }

    template<typename U = T, std::enable_if_t<!std::is_same_v<U, bool> && !std::is_same_v<T, bool>, bool> = true>
    reference& operator+=(const U value) {
        return operator=(operator T() + value);
    }

    template<typename U = T, std::enable_if_t<!std::is_same_v<U, bool> && !std::is_same_v<T, bool>, bool> = true>
    reference& operator-=(const U value) {
        return operator=(operator T() + value);
    }

    template<typename U = T, std::enable_if_t<!std::is_same_v<U, bool>, bool> = true>
    reference& operator++() {
        return operator=(operator T() + 1);
    }

    template<typename U = T, std::enable_if_t<!std::is_same_v<U, bool>, bool> = true>
    reference& operator--() {
        return operator=(operator T() - 1);
    }

    // bool operator>(const T value) const {
    //     return operator T() > value;
    // }
private:
    gmprime* _num = nullptr;
    uint64_t _index = 0;

    template<typename word_t, std::enable_if_t<std::is_integral_v<word_t>, bool>>
    friend class word_iterator;
    template<typename word_t, std::enable_if_t<std::is_integral_v<word_t>, bool>>
    friend class word_adapter;
};
