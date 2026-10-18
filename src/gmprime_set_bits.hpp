#pragma once

// Forward iterator over the positions of the bits set in a gmprime, in increasing order
class set_bit_iterator {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = size_t;
    using difference_type = std::ptrdiff_t;
    using pointer = const size_t*;
    using reference = const size_t&;

    set_bit_iterator() = default;

    set_bit_iterator(const gmprime* const num, const size_t pos): _num(num), _pos(pos) {
    }

    reference operator*() const {
        return _pos;
    }

    pointer operator->() const {
        return &_pos;
    }

    set_bit_iterator& operator++() {
        _pos = _num->next_set_bit(_pos + 1);
        return *this;
    }

    set_bit_iterator operator++(int) {
        const auto old = *this;
        ++*this;
        return old;
    }

    friend bool operator==(const set_bit_iterator& lhs, const set_bit_iterator& rhs) {
        return lhs._pos == rhs._pos;
    }

    friend bool operator!=(const set_bit_iterator& lhs, const set_bit_iterator& rhs) {
        return lhs._pos != rhs._pos;
    }

private:
    const gmprime* _num = nullptr;
    size_t _pos = gmprime::npos;
};

// Range of the positions of the bits set in a gmprime, see gmprime::set_bits
class set_bit_range {
public:
    explicit set_bit_range(const gmprime& num): _num(&num) {
    }

    // Keeps a temporary number alive as long as the range
    explicit set_bit_range(gmprime&& num): _owned(std::move(num)), _num(&*_owned) {
    }

    // Copies and moves point to their own copy of an owned number
    set_bit_range(const set_bit_range& other): _owned(other._owned), _num(_owned ? &*_owned : other._num) {
    }

    set_bit_range(set_bit_range&& other) noexcept: _owned(std::move(other._owned)), _num(_owned ? &*_owned : other._num) {
    }

    set_bit_range& operator=(const set_bit_range&) = delete;
    set_bit_range& operator=(set_bit_range&&) = delete;

    [[nodiscard]] set_bit_iterator begin() const {
        return {_num, _num->next_set_bit(0)};
    }

    [[nodiscard]] set_bit_iterator end() const {
        return {_num, gmprime::npos};
    }

private:
    std::optional<gmprime> _owned;
    const gmprime* _num;
};
