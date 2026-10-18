#ifndef GMPRIME_CODING_HPP
#define GMPRIME_CODING_HPP

// Online modeling and coding tools: predict the bits of a data stream from what precedes them, and code
// them with an arithmetic coder whose output is a gmprime. Every model is causal, so that a decoder
// rebuilds the same predictions from the bits it has already decoded.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>
#include "gmprime.hpp"

namespace GMPrime::coding {

// Logistic domain: stretch(p) = ln(p / (1 - p)), its inverse squash
inline double stretch(const double p) {
    return std::log(p / (1 - p));
}

inline double squash(const double x) {
    return 1 / (1 + std::exp(-x));
}

// Keeps a probability codable on 12 bits
inline double clamp_probability(const double p) {
    return std::min(std::max(p, 1.0 / 4096), 4095.0 / 4096);
}

// 64-bit hash finalizer (MurmurHash3): spreads close context values over a table
inline uint64_t hash_mix(uint64_t h) {
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdull;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ull;
    h ^= h >> 33;
    return h;
}

// Hash of the order bytes before position pos in history (missing ones count as 0), distinct per order
inline uint64_t context_hash(const std::vector<uint8_t>& history, const size_t pos, const size_t order) {
    uint64_t h = (order + 1) * 0x9E3779B97F4A7C15ull;
    for (size_t i = 1; i <= order; ++i) {
        h = (h ^ (pos >= i ? history[pos - i] : 0)) * 0x100000001B3ull;
    }
    return h;
}

// Adaptive probability that a bit is 1: moves toward each bit at rate 1 / (n + 1.5), fast while few bits
// were seen, then more slowly as n grows up to limit (a higher limit suits more stationary data)
struct bit_counter {
    float p = 0.5f;
    uint16_t n = 0;

    [[nodiscard]] bool seen() const {
        return n != 0;
    }

    void update(const int bit, const int limit) {
        p += static_cast<float>((bit - p) / (n + 1.5));
        if (n < limit) {
            ++n;
        }
    }
};

// Counts of the 0s and 1s seen in a context, in one byte: each bit increments its count (up to 15) and
// discounts the opposite one beyond 2, so that recent bits weigh more (nonstationary). The state stands
// for the history of the context, see history_map.
struct bit_history {
    uint8_t state = 0;

    [[nodiscard]] int zeros() const {
        return state >> 4;
    }

    [[nodiscard]] int ones() const {
        return state & 15;
    }

    void update(const int bit) {
        int same = bit ? ones() : zeros();
        int other = bit ? zeros() : ones();
        same = std::min(same + 1, 15);
        if (other > 2) {
            other = (other + 1) / 2;
        }
        const int n0 = bit ? other : same;
        const int n1 = bit ? same : other;
        state = static_cast<uint8_t>(n0 << 4 | n1);
    }
};

// Indirect modeling: probability of a 1 for each bit history state, learned across every context in that
// state. Rare contexts know little on their own; what usually follows "seen once, followed by a 1"
// generalizes from all of them.
class history_map {
public:
    explicit history_map(const int limit = 1023): _limit(limit) {
        for (int state = 0; state < 256; ++state) {
            // Starts from the counts themselves: (n1 + 1/2) / (n0 + n1 + 1)
            _counters[state].p = static_cast<float>(((state & 15) + 0.5) / ((state >> 4) + (state & 15) + 1.0));
        }
    }

    [[nodiscard]] double p(const bit_history history) const {
        return _counters[history.state].p;
    }

    void update(const bit_history history, const int bit) {
        _counters[history.state].update(bit, _limit);
    }

private:
    bit_counter _counters[256];
    int _limit;
};

// Slots of type T addressed by a context hash, 2^bits of them, without collision detection
template<typename T>
class hashed_table {
public:
    explicit hashed_table(const int bits): _slots(size_t{1} << bits), _mask((size_t{1} << bits) - 1) {
    }

    T& at(const uint64_t hash) {
        return _slots[hash_mix(hash) & _mask];
    }

private:
    std::vector<T> _slots;
    size_t _mask;
};

// Bit counters addressed by a context hash
using hashed_counters = hashed_table<bit_counter>;

// Everything a context knows about the next bit, in 8 bytes: a direct adaptive probability (as bit_counter),
// its bit history (for a history_map), and the checksum used by checked_table
struct context_slot {
    float p = 0.5f;
    uint16_t n = 0;
    bit_history history;
    uint8_t check = 0;

    [[nodiscard]] bool seen() const {
        return n != 0;
    }

    void update(const int bit, const int limit) {
        p += static_cast<float>((bit - p) / (n + 1.5));
        if (n < limit) {
            ++n;
        }
        history.update(bit);
    }
};

// The byte that followed the last occurrences of a context, and how many times in a row: predicts that it
// follows again (a run, or "what followed last time"). Has a check field for checked_table.
struct last_byte_slot {
    uint8_t byte = 0;
    uint8_t count = 0;
    uint8_t check = 0;

    void update(const uint8_t next) {
        if (count > 0 && byte == next) {
            count = static_cast<uint8_t>(std::min(count + 1, 255));
        } else {
            byte = next;
            count = 1;
        }
    }

    // Bit expected next given the bits of the current byte coded so far (c0 with a leading 1, bit_count of
    // them), or -1 if the remembered byte disagrees with them or nothing is remembered
    [[nodiscard]] int expected_bit(const uint32_t c0, const int bit_count) const {
        const int shift = 8 - bit_count;
        if (count == 0 || ((byte | 0x100) >> shift) != static_cast<int>(c0)) {
            return -1;
        }
        return (byte >> (shift - 1)) & 1;
    }
};

// Slots of type T (with a uint8_t check field) addressed by a context hash, with an 8-bit checksum of the
// hash in each slot: a context landing on a slot owned by another one resets it rather than sharing its
// statistics. The newest context wins. Fewer collisions than hashed_table, at 8 bits per slot.
template<typename T>
class checked_table {
public:
    explicit checked_table(const int bits): _slots(size_t{1} << bits), _mask((size_t{1} << bits) - 1) {
    }

    T& at(const uint64_t hash) {
        const uint64_t h = hash_mix(hash);
        auto& slot = _slots[h & _mask];
        // Never 0, so that a fresh slot always takes its first context's checksum
        const auto check = static_cast<uint8_t>((h >> 56) | 1);
        if (slot.check != check) {
            slot = T{};
            slot.check = check;
        }
        return slot;
    }

private:
    std::vector<T> _slots;
    size_t _mask;
};

// TEXT

// Class of a byte for selecting weights or contexts: what kind of character precedes the next one matters
// more for text than its exact value
enum class char_class : uint8_t { other, newline, space, lowercase, uppercase, digit, punctuation };

inline char_class classify(const uint8_t byte) {
    if (byte == '\n' || byte == '\r')
        return char_class::newline;
    if (byte == ' ' || byte == '\t')
        return char_class::space;
    if (byte >= 'a' && byte <= 'z')
        return char_class::lowercase;
    if (byte >= 'A' && byte <= 'Z')
        return char_class::uppercase;
    if (byte >= '0' && byte <= '9')
        return char_class::digit;
    if (byte >= 0x21 && byte <= 0x7E)
        return char_class::punctuation;
    return char_class::other;
}

// Line structure of a text, the 2D view of the bytes: the column of the next byte, the byte above it in the
// previous line, and the previous line's length (lines end with '\n')
class line_tracker {
public:
    // To call after appending a byte to history
    void byte_added(const std::vector<uint8_t>& history) {
        if (history.back() == '\n') {
            _previous_start = _start;
            _start = history.size();
        }
    }

    [[nodiscard]] size_t column(const std::vector<uint8_t>& history) const {
        return history.size() - _start;
    }

    [[nodiscard]] size_t previous_length() const {
        return _start - _previous_start;
    }

    // Byte at the same column in the previous line, 0 if that line is shorter or there is none
    [[nodiscard]] uint8_t above(const std::vector<uint8_t>& history) const {
        const size_t pos = _previous_start + column(history);
        return _start > 0 && pos < _start ? history[pos] : 0;
    }

private:
    size_t _start = 0, _previous_start = 0;
};

// Reversible capitalization transform: words are lowercased, a capitalized word (first letter uppercase, the
// others lowercase) preceded by 0x01, a word of several letters all uppercase by 0x02, other words kept as
// they are. "The" and "the" then share statistics. Returns nothing if text contains a 0x01 or 0x02 byte.
inline std::optional<std::vector<uint8_t>> capitals_encode(const std::vector<uint8_t>& text) {
    const auto upper = [](const uint8_t c) { return c >= 'A' && c <= 'Z'; };
    const auto lower = [](const uint8_t c) { return c >= 'a' && c <= 'z'; };
    std::vector<uint8_t> result;
    result.reserve(text.size() + text.size() / 16);

    for (size_t i = 0; i < text.size();) {
        if (text[i] == 0x01 || text[i] == 0x02) {
            return std::nullopt;
        }
        if (!upper(text[i]) && !lower(text[i])) {
            result.push_back(text[i++]);
            continue;
        }

        size_t end = i;
        bool all_upper = true, rest_lower = true;
        while (end < text.size() && (upper(text[end]) || lower(text[end]))) {
            all_upper = all_upper && upper(text[end]);
            rest_lower = rest_lower && (end == i || lower(text[end]));
            ++end;
        }

        if (all_upper && end - i >= 2) {
            result.push_back(0x02);
        } else if (upper(text[i]) && rest_lower) {
            result.push_back(0x01);
        } else {
            result.insert(result.end(), text.begin() + static_cast<std::ptrdiff_t>(i), text.begin() + static_cast<std::ptrdiff_t>(end));
            i = end;
            continue;
        }
        for (; i < end; ++i) {
            result.push_back(upper(text[i]) ? static_cast<uint8_t>(text[i] + ('a' - 'A')) : text[i]);
        }
    }

    return result;
}

inline std::vector<uint8_t> capitals_decode(const std::vector<uint8_t>& encoded) {
    // Every byte but the escapes comes back, recapitalized after one
    const auto escapes = std::count_if(encoded.begin(), encoded.end(), [](const uint8_t c) { return c == 0x01 || c == 0x02; });
    std::vector<uint8_t> result(encoded.size() - static_cast<size_t>(escapes));
    size_t out = 0;

    for (size_t i = 0; i < encoded.size(); ++i) {
        if (encoded[i] != 0x01 && encoded[i] != 0x02) {
            result[out++] = encoded[i];
            continue;
        }
        const bool whole_word = encoded[i] == 0x02;
        const size_t first = i + 1;
        while (i + 1 < encoded.size() && encoded[i + 1] >= 'a' && encoded[i + 1] <= 'z') {
            ++i;
            const bool capitalize = whole_word || i == first;
            result[out++] = capitalize ? static_cast<uint8_t>(encoded[i] - ('a' - 'A')) : encoded[i];
        }
    }

    return result;
}

// BIT GEOMETRY

// Bit p of bytes read as a stream, each byte from its highest bit: the order of a bitmap's pixels, where a
// fold of any width in bits lines up rows (gmprime's own bit order starts each byte from its lowest bit)
inline int stream_bit(const std::vector<uint8_t>& bytes, const size_t p) {
    return (bytes[p / 8] >> (7 - p % 8)) & 1;
}

// Packs the stream bits at pos - d for each distance d (at most 64, 0 before the start) into a context:
// with a fold of width L, distances 1, 2, L - 1, L, L + 1, 2L... make a 2D neighbourhood, L * H a 3D one
inline uint64_t bit_template(const std::vector<uint8_t>& bytes, const size_t pos, const std::vector<size_t>& distances) {
    uint64_t context = 0;
    for (const auto d : distances) {
        context = (context << 1) | static_cast<uint64_t>(d <= pos ? stream_bit(bytes, pos - d) : 0);
    }
    return context;
}

// Pixel runs at stream bit pos of a bitmap folded at width bits per row: left, the length of the run of equal
// bits ending just before pos (1 to limit); above, how many bits the one above pos stays the same toward the
// right (0 to limit). Bits before the start count as 0.
struct bitmap_runs {
    size_t left = 0;
    size_t above = 0;
};

inline bitmap_runs runs_at(const std::vector<uint8_t>& bytes, const size_t pos, const size_t width, const size_t limit = 64) {
    const auto bit = [&](const size_t distance) {
        return distance <= pos ? stream_bit(bytes, pos - distance) : 0;
    };
    bitmap_runs runs;
    runs.left = 1;
    while (runs.left < limit && bit(runs.left) == bit(runs.left + 1)) {
        ++runs.left;
    }
    while (runs.above < limit && runs.above + 2 < width && bit(width - runs.above) == bit(width - runs.above - 1)) {
        ++runs.above;
    }
    return runs;
}

// Size of a context table for data of the given size, in bits of slot count: room for about four slots per
// bit of data, between min_bits and max_bits (22, 32 MB of context_slot, suits data from 128 KB)
inline int table_bits(const size_t data_bytes, const int max_bits = 22, const int min_bits = 16) {
    int bits = 0;
    while (bits < max_bits && (size_t{1} << bits) < data_bytes * 8 * 4) {
        ++bits;
    }
    return std::max(bits, min_bits);
}

// Cost of each byte coded with a model, -log2 of the probabilities given to its actual bits: where the bits
// go, to compare models or study the data
class cost_profile {
public:
    // Records a coded bit and the probability that it is 1 given to the coder
    void add(const int bit, const double p) {
        _current -= std::log2(bit ? p : 1 - p);
        if (++_bits == 8) {
            _costs.push_back(_current);
            _current = 0;
            _bits = 0;
        }
    }

    // Bits spent on each complete byte
    [[nodiscard]] const std::vector<double>& per_byte() const {
        return _costs;
    }

    [[nodiscard]] double total() const {
        double sum = _current;
        for (const auto cost : _costs) {
            sum += cost;
        }
        return sum;
    }

private:
    std::vector<double> _costs;
    double _current = 0;
    int _bits = 0;
};

// Online logistic mixer: p = squash(sum of w_i x_i) over stretched input probabilities x_i, with weights
// learned by gradient descent on the coding cost, one weight vector per selected set (e.g. per context)
class logistic_mixer {
public:
    logistic_mixer(const size_t inputs, const size_t sets, const double rate, const double initial_weight = 0.3)
        : _weights(inputs * sets, initial_weight), _inputs(inputs), _rate(rate) {
        _x.reserve(inputs);
    }

    void add(const double x) {
        _x.push_back(x);
    }

    // Adds several inputs at once, in order: cheaper than one add per input when they are many
    void add(const std::vector<double>& inputs) {
        _x.insert(_x.end(), inputs.begin(), inputs.end());
    }

    // Mixes the inputs added since the last update with the weights of set
    double mix(const size_t set) {
        _set = set;
        double dot = 0;
        for (size_t i = 0; i < _x.size(); ++i) {
            dot += _weights[_set * _inputs + i] * _x[i];
        }
        _p = clamp_probability(squash(dot));
        return _p;
    }

    void update(const int bit) {
        const double err = bit - _p;
        for (size_t i = 0; i < _x.size(); ++i) {
            _weights[_set * _inputs + i] += _rate * err * _x[i];
        }
        _x.clear();
    }

    // Changes the learning rate, e.g. for a schedule decreasing with the bits seen
    void set_rate(const double rate) {
        _rate = rate;
    }

private:
    std::vector<double> _weights, _x;
    size_t _inputs, _set = 0;
    double _rate, _p = 0.5;
};

// Adaptive probability map (secondary estimation): refines a probability given a context, interpolating
// between 33 learned values over its stretched range
class probability_map {
public:
    probability_map(const size_t contexts, const double rate): _values(contexts * 33), _rate(rate) {
        for (size_t i = 0; i < _values.size(); ++i) {
            _values[i] = static_cast<float>(squash((static_cast<int>(i % 33) - 16) / 2.0));
        }
    }

    double refine(const double p, const size_t context) {
        const double s = std::max(-7.999, std::min(7.999, stretch(p))) * 2 + 16;
        const auto low = static_cast<size_t>(s);
        _weight = s - static_cast<double>(low);
        _index = context * 33 + low;
        return clamp_probability(_values[_index] * (1 - _weight) + _values[_index + 1] * _weight);
    }

    void update(const int bit) {
        _values[_index] += static_cast<float>((bit - _values[_index]) * _rate * (1 - _weight));
        _values[_index + 1] += static_cast<float>((bit - _values[_index + 1]) * _rate * _weight);
    }

private:
    std::vector<float> _values;
    double _rate, _weight = 0;
    size_t _index = 0;
};

// Causal match model: after each byte, looks for an earlier occurrence of the last min_length bytes and
// predicts that the bytes which followed it come again, bit by bit, until a prediction fails
class match_model {
public:
    explicit match_model(const size_t min_length, const int table_bits = 22)
        : _table(size_t{1} << table_bits, 0), _min_length(min_length) {
    }

    // To call after appending a byte to history
    void byte_added(const std::vector<uint8_t>& history) {
        const size_t pos = history.size();
        if (_length > 0 && _ptr + 1 < pos) {
            ++_ptr;
            ++_length;
        } else {
            _length = 0;
        }
        if (pos < _min_length) {
            return;
        }

        uint64_t h = 0;
        for (size_t i = 1; i <= _min_length; ++i) {
            h = (h ^ history[pos - i]) * 0x100000001B3ull;
        }
        auto& slot = _table[hash_mix(h) & (_table.size() - 1)];

        if (_length == 0 && slot > 0) {
            // Checks the candidate, the table holding no key, and measures how far back it matches
            size_t length = 0;
            while (length < slot && length < 65535 && history[slot - 1 - length] == history[pos - 1 - length]) {
                ++length;
            }
            if (length >= _min_length) {
                _ptr = slot;
                _length = static_cast<uint32_t>(length);
            }
        }
        slot = static_cast<uint32_t>(pos);
    }

    // Bit expected next, given the bits of the current byte coded so far (c0 with a leading 1, bit_count
    // of them), or -1 without a prediction
    [[nodiscard]] int expected_bit(const std::vector<uint8_t>& history, const uint32_t c0, const int bit_count) const {
        if (_length == 0) {
            return -1;
        }
        const int predicted = history[_ptr];
        const int shift = 8 - bit_count;
        if (((predicted | 0x100) >> shift) != static_cast<int>(c0)) {
            return -1;
        }
        return (predicted >> (shift - 1)) & 1;
    }

    // Length of the current match in bytes, 0 if none
    [[nodiscard]] uint32_t length() const {
        return _length;
    }

    // To call when a predicted bit turned out wrong: the match ends
    void miss() {
        _length = 0;
    }

private:
    std::vector<uint32_t> _table;
    size_t _min_length, _ptr = 0;
    uint32_t _length = 0;
};

// Binary arithmetic coder on 32 bits (carryless), probabilities rounded to 12 bits. The code is the gmprime
// whose bytes are the output, the first one lowest: the coded message as one number.
inline uint32_t probability_12(const double p) {
    return static_cast<uint32_t>(std::min(4095.0, std::max(1.0, std::round(p * 4096))));
}

class arithmetic_encoder {
public:
    // Codes bit, p being the probability that it is 1
    void encode(const int bit, const double p) {
        const uint32_t mid = _x1 + static_cast<uint32_t>((static_cast<uint64_t>(_x2 - _x1) * probability_12(p)) >> 12);
        if (bit) {
            _x2 = mid;
        } else {
            _x1 = mid + 1;
        }
        while (((_x1 ^ _x2) & 0xFF000000) == 0) {
            _out.push_back(static_cast<uint8_t>(_x2 >> 24));
            _x1 <<= 8;
            _x2 = (_x2 << 8) | 255;
        }
    }

    [[nodiscard]] gmprime finish() {
        for (int i = 0; i < 4; ++i) {
            _out.push_back(static_cast<uint8_t>(_x1 >> 24));
            _x1 <<= 8;
        }
        return gmprime::from_bytes(_out);
    }

private:
    std::vector<uint8_t> _out;
    uint32_t _x1 = 0, _x2 = 0xFFFFFFFF;
};

class arithmetic_decoder {
public:
    explicit arithmetic_decoder(const gmprime& code): _in(code.bytes()) {
        for (int i = 0; i < 4; ++i) {
            _x = (_x << 8) | next();
        }
    }

    // Decodes a bit, p being the probability that it is 1, as given to the encoder
    int decode(const double p) {
        const uint32_t mid = _x1 + static_cast<uint32_t>((static_cast<uint64_t>(_x2 - _x1) * probability_12(p)) >> 12);
        const int bit = _x <= mid ? 1 : 0;
        if (bit) {
            _x2 = mid;
        } else {
            _x1 = mid + 1;
        }
        while (((_x1 ^ _x2) & 0xFF000000) == 0) {
            _x1 <<= 8;
            _x2 = (_x2 << 8) | 255;
            _x = (_x << 8) | next();
        }
        return bit;
    }

private:
    uint8_t next() {
        return _pos < _in.size() ? _in[_pos++] : 0;
    }

    std::vector<unsigned char> _in;
    size_t _pos = 0;
    uint32_t _x1 = 0, _x2 = 0xFFFFFFFF, _x = 0;
};

} // GMPrime::coding

#endif //GMPRIME_CODING_HPP
