#pragma once

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t>, bool> = true>
class word_iterator {
public:
    using iterator_category = std::bidirectional_iterator_tag;

    word_iterator() = default;
    word_iterator(const word_iterator& other) = default;
    word_iterator(word_iterator&& other) = default;

    explicit word_iterator(gmprime& value, const int64_t index, const bool reverse = false) {
        _data = value.operator[]<word_t>(std::max(static_cast<int64_t>(0), index));
        _index = index;
        _reverse = reverse;
    }

    word_iterator& operator=(const word_iterator& rhs) {
        if (this != &rhs) {
            _data = rhs._data;
            _index = rhs._index;
            _reverse = rhs._reverse;
        }
        return *this;
    }

    constexpr reference<word_t>& operator*() noexcept {
        return *_current;
    }

    constexpr reference<word_t>& operator*() const noexcept {
        return *_current;
    }

    reference<word_t>* operator->() {
        return _current;
    }

    reference<word_t>* operator->() const {
        return _current;
    }

    word_iterator& operator++() {
        if (_reverse) {
            --_current->_index;
            --_index;
        } else {
            ++_current->_index;
            ++_index;
        }
        return *this;
    }

    word_iterator& operator--() {
        if (_reverse) {
            ++_current->_index;
            ++_index;
        } else {
            --_current->_index;
            --_index;
        }
        return *this;
    }

    word_iterator operator++(int) {
        const word_iterator tmp = *this;
        ++*this;
        return tmp;
    }

    word_iterator operator--(int) {
        const word_iterator tmp = *this;
        --*this;
        return tmp;
    }

    bool operator==(const word_iterator& other) const {
        return _index == other._index && _current->_num == other._current->_num;
    }

    bool operator!=(const word_iterator& other) const {
        return _index != other._index || _current->_num != other._current->_num;
    }
private:
    reference<word_t> _data;
    reference<word_t>* const _current = &_data;
    int64_t _index = 0;
    bool _reverse = false;
};

template<typename word_t, std::enable_if_t<std::is_integral_v<word_t>, bool>>
class word_adapter {
public:
    using iterator = word_iterator<word_t>;

    explicit word_adapter(gmprime& value, const bool reverse = false) : _value(value), _reverse(reverse) {
    }

    [[nodiscard]] constexpr iterator begin() const noexcept {
        const auto index = _reverse ? size() - 1 : 0;
        return iterator(_value, index, _reverse);
    }

    [[nodiscard]] constexpr iterator end() const noexcept {
        const int64_t index = _reverse ? -1 : size();
        return iterator(_value, index, _reverse);
    }

    [[nodiscard]] constexpr size_t size() const noexcept {
        if constexpr (std::is_same_v<word_t, bool>) {
            return _value.size_in_base(2);
        } else {
            return _value.used_words<word_t>();
        }
    }
private:
    gmprime& _value;
    bool _reverse = false;
};
