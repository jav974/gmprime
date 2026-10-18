#ifndef GMPRIME_CODING_HPP
#define GMPRIME_CODING_HPP

// Online modeling and coding tools: predict the bits of a data stream from what precedes them, and code
// them with an arithmetic coder whose output is a gmprime. Every model is causal, so that a decoder
// rebuilds the same predictions from the bits it has already decoded.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <new>
#include <optional>
#include <type_traits>
#include <vector>
#include "gmprime.hpp"

namespace GMPrime::coding {

// INTEGER ARITHMETIC: probabilities on 16 bits (0 to 65535 for 0 to 1; counters keep 22 bits), the logistic
// domain in 1/256 units over [-12, 12] (-3072 to 3072), mixer weights in 8.24 fixed point, rates scaled by 65536.
// No floating point in the models: tables come from integer code or constants, so every platform and compiler
// computes exactly the same predictions, and a message decodes everywhere (floating point does not guarantee
// it: compilers may fuse multiplications and additions differently).

// Products of negative values are shifted right throughout, which must keep their sign (guaranteed from C++20,
// done by every mainstream compiler before)
static_assert((-5 >> 1) == -3 && (int64_t{-5} >> 1) == -3, "right shifts must be arithmetic");

// squash(d / 256) * 65536 at d = -3072, -3040, ..., 3072, rounded (constants: no exp at run time)
inline const uint16_t* squash_points() {
    static const uint16_t points[193] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
    2, 2, 2, 3, 3, 3, 4, 4, 5, 6, 6, 7,
    8, 9, 10, 12, 13, 15, 17, 19, 22, 25, 28, 32,
    36, 41, 47, 53, 60, 68, 77, 87, 98, 111, 126, 143,
    162, 184, 208, 236, 267, 302, 342, 387, 439, 497, 562, 636,
    720, 815, 922, 1042, 1179, 1333, 1506, 1701, 1921, 2168, 2446, 2758,
    3108, 3500, 3938, 4427, 4971, 5577, 6249, 6992, 7812, 8714, 9702, 10782,
    11955, 13226, 14595, 16062, 17625, 19282, 21025, 22849, 24743, 26695, 28693, 30723,
    32768, 34813, 36843, 38841, 40793, 42687, 44511, 46254, 47911, 49474, 50941, 52310,
    53581, 54754, 55834, 56822, 57724, 58544, 59287, 59959, 60565, 61109, 61598, 62036,
    62428, 62778, 63090, 63368, 63615, 63835, 64030, 64203, 64357, 64494, 64614, 64721,
    64816, 64900, 64974, 65039, 65097, 65149, 65194, 65234, 65269, 65300, 65328, 65352,
    65374, 65393, 65410, 65425, 65438, 65449, 65459, 65468, 65476, 65483, 65489, 65495,
    65500, 65504, 65508, 65511, 65514, 65517, 65519, 65521, 65523, 65524, 65526, 65527,
    65528, 65529, 65530, 65530, 65531, 65532, 65532, 65533, 65533, 65533, 65534, 65534,
    65534, 65534, 65535, 65535, 65535, 65535, 65535, 65535, 65535, 65535, 65535, 65535,
    65535,
    };
    return points;
}

// 1 / (1 + e^(-d / 256)) on 16 bits, linearly interpolated between the constant points, clamped beyond 12
inline int squash(int d) {
    if (d >= 3072) return 65535;
    if (d <= -3072) return 1;
    const int i = (d + 3072) >> 5;
    const int w = (d + 3072) & 31;
    const uint16_t* p = squash_points();
    return (p[i] * (32 - w) + p[i + 1] * w + 16) >> 5;
}

// ln(p / (1 - p)) * 256 for p on 16 bits: the inverse of squash, as a table built from it
inline int stretch(const int p16) {
    static const auto table = [] {
        std::vector<int16_t> t(65536);
        int next = 0;
        for (int d = -3071; d <= 3071; ++d) {
            const int v = squash(d);
            for (int i = next; i <= v && i < 65536; ++i) {
                t[static_cast<size_t>(i)] = static_cast<int16_t>(d);
            }
            next = std::max(next, v + 1);
        }
        for (int i = next; i < 65536; ++i) {
            t[static_cast<size_t>(i)] = 3071;
        }
        return t;
    }();
    return table[static_cast<size_t>(std::min(65535, std::max(0, p16)))];
}

// Keeps a 16-bit probability within the coder's 12 bits: 1/4096 to 4095/4096
inline int clamp_probability(const int p16) {
    return std::min(65520, std::max(16, p16));
}

// The coder's 12-bit probability of a 16-bit one
inline uint32_t probability_12(const int p16) {
    return static_cast<uint32_t>(std::min(4095, std::max(1, (p16 + 8) >> 4)));
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

// 65536 / (n + 1.5) for n up to 1023, the learning rates of the counters, by integer division
inline int counter_rate(const int n) {
    static const auto table = [] {
        std::vector<int> values(1024);
        for (int i = 0; i < 1024; ++i) {
            values[static_cast<size_t>(i)] = 131072 / (2 * i + 3);
        }
        return values;
    }();
    return table[static_cast<size_t>(n)];
}

// Moves a 22-bit probability toward a bit at rate counter_rate(n) / 65536. 22 bits, not 16: with slow rates
// (n large), steps on 16 bits round to 0 well before the probability nears what the coder can use
inline void learn(uint32_t& p, const int bit, const int n) {
    const int64_t target = bit ? (1 << 22) - 1 : 0;
    p = static_cast<uint32_t>(static_cast<int64_t>(p) + (((target - static_cast<int64_t>(p)) * counter_rate(n)) >> 16));
}

// Adaptive probability that a bit is 1, on 16 bits: moves toward each bit at rate 1 / (n + 1.5), fast while
// few bits were seen, then more slowly as n grows up to limit
struct bit_counter {
    uint32_t p = 1 << 21;   // on 22 bits
    uint16_t n = 0;

    [[nodiscard]] bool seen() const {
        return n != 0;
    }

    // The probability on 16 bits
    [[nodiscard]] int p16() const {
        return static_cast<int>(p >> 6);
    }

    void update(const int bit, const int limit) {
        learn(p, bit, n);
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
            const int n0 = state >> 4, n1 = state & 15;
            _counters[state].p = static_cast<uint32_t>((2 * n1 + 1) * (1 << 22) / (2 * (n0 + n1 + 1)));
        }
    }

    [[nodiscard]] int p(const bit_history history) const {
        return _counters[history.state].p16();
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
    uint32_t p = 1 << 21;   // on 22 bits
    uint16_t n = 0;
    bit_history history;
    uint8_t check = 0;

    [[nodiscard]] bool seen() const {
        return n != 0;
    }

    [[nodiscard]] int p16() const {
        return static_cast<int>(p >> 6);
    }

    void update(const int bit, const int limit) {
        learn(p, bit, n);
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
    // Zeroed memory from calloc: a slot of zero bytes holds checksum 0, which no context has, so it reads as
    // empty, and the system provides the pages only once touched (fast to create, light on small inputs)
    explicit checked_table(const int bits)
        : _slots(static_cast<T*>(std::calloc(size_t{1} << bits, sizeof(T)))), _mask((size_t{1} << bits) - 1) {
        static_assert(std::is_trivially_copyable_v<T>, "slots are created from zeroed memory");
        if (!_slots) {
            throw std::bad_alloc(); // LCOV_EXCL_LINE: sanitizers abort on such a request rather than return null
        }
    }

    T& at(const uint64_t hash) {
        return at(hash, 0, ~size_t{0});
    }

    // Asks the processor to bring the slot of hash into cache, without reading it: to overlap the memory
    // latency of the next accesses with the work on the current ones
    void prefetch(const uint64_t hash) const {
#if defined(__GNUC__) || defined(__clang__)
        __builtin_prefetch(&_slots.get()[hash_mix(hash) & _mask]);
#else
        static_cast<void>(hash);
#endif
    }

    // Slot offset (below 16) of a bucket of 16 adjacent slots shared by one hash: e.g. the nodes of a nibble's
    // bit tree, so that coding a nibble touches one bucket instead of a slot anywhere for each bit
    T& at(const uint64_t hash, const size_t offset) {
        return at(hash, offset, ~size_t{15});
    }

private:
    T& at(const uint64_t hash, const size_t offset, const size_t base_mask) {
        const uint64_t h = hash_mix(hash);
        auto& slot = _slots.get()[(h & _mask & base_mask) | offset];
        // Never 0, so that a fresh slot always takes its first context's checksum
        const auto check = static_cast<uint8_t>((h >> 56) | 1);
        if (slot.check != check) {
            slot = T{};
            slot.check = check;
        }
        return slot;
    }

    struct free_deleter {
        void operator()(T* slots) const {
            std::free(slots);
        }
    };

    std::unique_ptr<T, free_deleter> _slots;
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

// Stream bits [first, first + count), count from 1 to 57, packed the first one highest: a word at a time instead
// of bit by bit. first may be negative: bits before the start, and beyond the end, are 0
inline uint64_t stream_bits(const std::vector<uint8_t>& bytes, const int64_t first, const int count) {
    if (first + count <= 0) {
        return 0;
    }
    if (first < 0) {
        return stream_bits(bytes, 0, static_cast<int>(first + count));
    }
    const auto byte = static_cast<size_t>(first) / 8;
    const auto offset = static_cast<int>(first % 8);
    uint64_t word = 0;
    for (size_t i = 0; i < 8; ++i) {
        word = (word << 8) | (byte + i < bytes.size() ? bytes[byte + i] : 0);
    }
    return (word << offset) >> (64 - count);
}

// The low count bits of value in reverse order
inline uint64_t reverse_bits(uint64_t value, const int count) {
    uint64_t reversed = 0;
    for (int i = 0; i < count; ++i) {
        reversed = (reversed << 1) | (value & 1);
        value >>= 1;
    }
    return reversed;
}

// bit_template with its distances grouped once into runs of consecutive ones (d, d + 1... or d, d - 1...), each
// read as one word: the same context in a few word reads instead of one bit read per distance, for templates
// used at every bit
class compiled_template {
public:
    explicit compiled_template(const std::vector<size_t>& distances) {
        for (size_t i = 0; i < distances.size();) {
            segment s{distances[i], 1, false};
            if (i + 1 < distances.size() && (distances[i + 1] == distances[i] + 1 || distances[i + 1] + 1 == distances[i])) {
                s.ascending = distances[i + 1] == distances[i] + 1;
                while (i + static_cast<size_t>(s.length) < distances.size() && s.length < 57
                       && distances[i + static_cast<size_t>(s.length)] == (s.ascending ? distances[i] + static_cast<size_t>(s.length) : distances[i] - static_cast<size_t>(s.length))) {
                    ++s.length;
                }
            }
            _segments.push_back(s);
            i += static_cast<size_t>(s.length);
        }
    }

    [[nodiscard]] uint64_t at(const std::vector<uint8_t>& bytes, const size_t pos) const {
        uint64_t context = 0;
        for (const auto& s : _segments) {
            const auto p = static_cast<int64_t>(pos) - static_cast<int64_t>(s.distance);
            // Descending distances are ascending positions, in stream order; ascending ones the reverse
            const uint64_t bits = s.ascending ? reverse_bits(stream_bits(bytes, p - s.length + 1, s.length), s.length)
                                              : stream_bits(bytes, p, s.length);
            context = (context << s.length) | bits;
        }
        return context;
    }

private:
    struct segment {
        size_t distance;
        int length;
        bool ascending;
    };
    std::vector<segment> _segments;
};

// Number of zero bits below the lowest one, and above the highest one, of a non-zero word
inline int trailing_zeros(uint64_t word) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(word);
#else
    int count = 0;
    for (; (word & 1) == 0; word >>= 1) {
        ++count;
    }
    return count;
#endif
}

inline int leading_zeros(uint64_t word) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(word);
#else
    int count = 0;
    for (; (word >> 63) == 0; word <<= 1) {
        ++count;
    }
    return count;
#endif
}

// Pixel runs at stream bit pos of a bitmap folded at width bits per row: left, the length of the run of equal
// bits ending just before pos (1 to limit); above, how many bits the one above pos stays the same toward the
// right (0 to limit). Bits before the start count as 0.
struct bitmap_runs {
    size_t left = 0;
    size_t above = 0;
};

inline bitmap_runs runs_at(const std::vector<uint8_t>& bytes, const size_t pos, const size_t width, const size_t limit = 64) {
    // Bits compared 56 pairs at a time: in a word of 57 stream bits, word ^ (word >> 1) marks where neighbours differ
    constexpr uint64_t pairs = (uint64_t{1} << 56) - 1;
    const auto p = static_cast<int64_t>(pos);
    bitmap_runs runs;

    // Left: equal neighbours going back from pos - 1, the lowest bit of each word
    runs.left = 1;
    for (int64_t end = p - 1; runs.left < limit; end -= 56) {
        const uint64_t word = stream_bits(bytes, end - 56, 57);
        const uint64_t differ = (word ^ (word >> 1)) & pairs;
        const auto equal = static_cast<size_t>(differ == 0 ? 56 : trailing_zeros(differ));
        runs.left = std::min(limit, runs.left + equal);
        if (equal < 56) {
            break;
        }
    }

    // Above: equal neighbours going forward from pos - width, the highest bit of each word, before pos - 1
    const size_t cap = std::min(limit, width >= 2 ? width - 2 : 0);
    for (int64_t first = p - static_cast<int64_t>(width); runs.above < cap; first += 56) {
        const uint64_t word = stream_bits(bytes, first, 57);
        const uint64_t differ = (word ^ (word >> 1)) & pairs;
        const auto equal = static_cast<size_t>(differ == 0 ? 56 : leading_zeros(differ) - 8);
        runs.above = std::min(cap, runs.above + equal);
        if (equal < 56) {
            break;
        }
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

// CRC-32 of bytes (the one of zip and png), to check that a message decoded to the original data. crc continues
// a previous checksum: crc32(b, crc32(a)) is the checksum of a followed by b
inline uint32_t crc32(const std::vector<uint8_t>& bytes, uint32_t crc = 0) {
    static const auto table = [] {
        std::vector<uint32_t> values(256);
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int k = 0; k < 8; ++k) {
                c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            }
            values[i] = c;
        }
        return values;
    }();
    crc = ~crc;
    for (const auto byte : bytes) {
        crc = table[(crc ^ byte) & 0xFF] ^ (crc >> 8);
    }
    return ~crc;
}

// Cost of each byte coded with a model, -log2 of the probabilities given to its actual bits: where the bits
// go, to compare models or study the data
class cost_profile {
public:
    // Records a coded bit and the probability that it is 1 given to the coder
    void add(const int bit, const int p16) {
        const double p = p16 / 65536.0;
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
    // Weights in 8.24 fixed point (24 bits of fraction, so that the tiny errors of near-certain predictions
    // still move them). rate: the learning rate times 65536, e.g. 197 for 0.003; weights start at
    // initial_weight / 65536
    logistic_mixer(const size_t inputs, const size_t sets, const int rate, const int initial_weight = 19661)
        : _weights(inputs * sets, initial_weight * 256), _inputs(inputs), _rate(rate) {
        _x.reserve(inputs);
    }

    // An input in the logistic domain, in 1/256 units (a stretched probability)
    void add(const int x) {
        _x.push_back(x);
    }

    // Mixes the inputs added since the last update with the weights of set: a 16-bit probability
    int mix(const size_t set) {
        const int p = mix(_x, set);
        _in = nullptr;
        return p;
    }

    // Mixes inputs kept by the caller, without copying them: they must stay unchanged until update
    int mix(const std::vector<int>& inputs, const size_t set) {
        _set = set;
        _in = &inputs;
        const int32_t* w = &_weights[_set * _inputs];
        int64_t dot = 0;
        for (size_t i = 0; i < inputs.size(); ++i) {
            dot += static_cast<int64_t>(w[i]) * inputs[i];
        }
        _p = clamp_probability(squash(static_cast<int>(std::max<int64_t>(-3072, std::min<int64_t>(3072, dot >> 24)))));
        return _p;
    }

    void update(const int bit) {
        const int64_t err = (bit ? 65535 : 0) - _p;
        int32_t* w = &_weights[_set * _inputs];
        const auto& x = _in ? *_in : _x;
        for (size_t i = 0; i < x.size(); ++i) {
            w[i] += static_cast<int32_t>((x[i] * err * _rate) >> 16);
        }
        _x.clear();
    }

    void set_rate(const int rate) {
        _rate = rate;
    }

private:
    std::vector<int32_t> _weights;
    std::vector<int> _x;
    const std::vector<int>* _in = nullptr;
    size_t _inputs, _set = 0;
    int64_t _rate;
    int _p = 32768;
};

// Adaptive probability map (secondary estimation): refines a 16-bit probability given a context, interpolating
// between 33 learned values over its stretched range [-8, 8]; rate in 1/65536 (e.g. 1311 for 0.02)
class probability_map {
public:
    probability_map(const size_t contexts, const int rate): _values(contexts * 33), _rate(rate) {
        for (size_t i = 0; i < 33 && i < _values.size(); ++i) {
            _values[i] = static_cast<uint32_t>(squash((static_cast<int>(i) - 16) * 128)) << 6;
        }
        for (size_t i = 33; i < _values.size(); ++i) {
            _values[i] = _values[i - 33];
        }
    }

    int refine(const int p16, const size_t context) {
        const int s = std::max(-2047, std::min(2047, stretch(p16))) + 2048;
        _weight = s & 127;
        _index = context * 33 + static_cast<size_t>(s >> 7);
        return clamp_probability(static_cast<int>(((static_cast<int64_t>(_values[_index]) * (128 - _weight) + static_cast<int64_t>(_values[_index + 1]) * _weight) >> 7) >> 6));
    }

    void update(const int bit) {
        const int64_t target = bit ? (1 << 22) - 1 : 0;
        _values[_index] = static_cast<uint32_t>(_values[_index] + (((target - _values[_index]) * _rate * (128 - _weight)) >> 23));
        _values[_index + 1] = static_cast<uint32_t>(_values[_index + 1] + (((target - _values[_index + 1]) * _rate * _weight) >> 23));
    }

private:
    std::vector<uint32_t> _values;   // on 22 bits
    int64_t _rate;
    int _weight = 0;
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

class arithmetic_encoder {
public:
    // Codes bit, p being the probability that it is 1
    void encode(const int bit, const int p) {
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
    int decode(const int p) {
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
