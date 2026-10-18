#include <cmath>
#include <unordered_map>
#include "gmprime.hpp"

using namespace GMPrime;

namespace {

// Hash table counting sequences of symbols, keyed by their start in a symbol array and compared by
// content: no copy of each context
class sequence_counts {
public:
    sequence_counts(const std::vector<uint32_t> &symbols, const size_t length)
        : _counts(0, hasher{&symbols, length}, equal{&symbols, length}) {
    }

    size_t &operator[](const size_t start) {
        return _counts[start];
    }

    [[nodiscard]] size_t at(const size_t start) const {
        return _counts.at(start);
    }

    [[nodiscard]] size_t size() const {
        return _counts.size();
    }

    [[nodiscard]] auto begin() const {
        return _counts.begin();
    }

    [[nodiscard]] auto end() const {
        return _counts.end();
    }

private:
    struct hasher {
        const std::vector<uint32_t> *symbols;
        size_t length;

        // FNV-1a over the symbols
        size_t operator()(const size_t start) const {
            uint64_t hash = 14695981039346656037ull;
            for (size_t i = start; i < start + length; ++i) {
                hash = (hash ^ (*symbols)[i]) * 1099511628211ull;
            }
            return static_cast<size_t>(hash ^ (hash >> 32));
        }
    };

    struct equal {
        const std::vector<uint32_t> *symbols;
        size_t length;

        bool operator()(const size_t a, const size_t b) const {
            const auto first = symbols->begin();
            return std::equal(first + a, first + a + length, first + b);
        }
    };

    std::unordered_map<size_t, size_t, hasher, equal> _counts;
};

void check_symbol_width(const size_t width) {
    if (width == 0 || width > 32) {
        throw std::invalid_argument("Symbol width must be between 1 and 32");
    }
}

// Suffix array by prefix doubling with radix sorts: O(n log n), ranks of the prefixes of length 2k
// computed from those of length k
std::vector<uint32_t> suffix_array(const std::vector<unsigned char> &text) {
    const size_t n = text.size();
    std::vector<uint32_t> sa(n), rank(n), next_rank(n), by_second(n);

    if (n == 0) {
        return sa;
    }

    std::vector<size_t> count(257, 0);
    for (const auto c : text) {
        ++count[c + 1];
    }
    for (size_t c = 1; c < count.size(); ++c) {
        count[c] += count[c - 1];
    }
    for (size_t i = 0; i < n; ++i) {
        sa[count[text[i]]++] = static_cast<uint32_t>(i);
    }

    rank[sa[0]] = 0;
    for (size_t r = 1; r < n; ++r) {
        rank[sa[r]] = rank[sa[r - 1]] + (text[sa[r]] != text[sa[r - 1]] ? 1 : 0);
    }

    for (size_t k = 1; rank[sa[n - 1]] < n - 1; k *= 2) {
        // Sorted by second half: suffixes without one first (shorter than k, so all distinct), then the others
        size_t p = 0;
        for (size_t i = n > k ? n - k : 0; i < n; ++i) {
            by_second[p++] = static_cast<uint32_t>(i);
        }
        for (size_t r = 0; r < n; ++r) {
            if (sa[r] >= k) {
                by_second[p++] = static_cast<uint32_t>(sa[r] - k);
            }
        }

        // Then stably by first half
        count.assign(static_cast<size_t>(rank[sa[n - 1]]) + 2, 0);
        for (size_t i = 0; i < n; ++i) {
            ++count[rank[i] + 1];
        }
        for (size_t c = 1; c < count.size(); ++c) {
            count[c] += count[c - 1];
        }
        for (size_t i = 0; i < n; ++i) {
            sa[count[rank[by_second[i]]]++] = by_second[i];
        }

        // Equal ranks for equal (first half, second half) pairs, a missing second half being distinct
        const auto second = [&](const size_t i) -> int64_t { return i + k < n ? rank[i + k] : -1; };
        next_rank[sa[0]] = 0;
        for (size_t r = 1; r < n; ++r) {
            const bool same = rank[sa[r]] == rank[sa[r - 1]] && second(sa[r]) == second(sa[r - 1]);
            next_rank[sa[r]] = next_rank[sa[r - 1]] + (same ? 0 : 1);
        }
        rank.swap(next_rank);
    }

    return sa;
}

// lcp[r]: length of the common prefix of the suffixes sa[r - 1] and sa[r], 0 for r = 0 (Kasai et al., O(n))
std::vector<uint32_t> lcp_array(const std::vector<unsigned char> &text, const std::vector<uint32_t> &sa) {
    const size_t n = text.size();
    std::vector<uint32_t> rank(n), lcp(n, 0);

    for (size_t r = 0; r < n; ++r) {
        rank[sa[r]] = static_cast<uint32_t>(r);
    }

    // The common prefix with the previous suffix in order shrinks by at most one from i to i + 1
    size_t h = 0;
    for (size_t i = 0; i < n; ++i) {
        if (rank[i] == 0) {
            h = 0;
            continue;
        }
        const size_t j = sa[rank[i] - 1];
        while (i + h < n && j + h < n && text[i + h] == text[j + h]) {
            ++h;
        }
        lcp[rank[i]] = static_cast<uint32_t>(h);
        if (h > 0) {
            --h;
        }
    }

    return lcp;
}

} // namespace

size_t gmprime::analyzed_bits() const {
    // Whole bytes, like a file: the high byte of a text is partial (ASCII leaves its top bit clear)
    const auto used_bytes = (used_bits() + CHAR_BIT - 1) / CHAR_BIT;
    return std::max(used_bytes, _byte_size) * CHAR_BIT;
}

std::vector<uint32_t> gmprime::symbols(const size_t width) const {
    check_symbol_width(width);

    std::vector<uint32_t> result(analyzed_bits() / width);
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = read_bits<uint32_t>(i * width, width);
    }

    return result;
}

gmprime gmprime::from_bytes(const std::vector<unsigned char> &bytes) {
    gmprime result;
    // GMP's C functions never throw (it aborts on allocation failure): no cleanup to plan for
    [&]() noexcept {
        if (!bytes.empty()) {
            mpz_import(result._value.get_mpz_t(), bytes.size(), -1, 1, 0, 0, bytes.data());
        }
    }();
    result._byte_size = bytes.size();
    return result;
}

std::vector<unsigned char> gmprime::bytes() const {
    std::vector<unsigned char> result(analyzed_bits() / CHAR_BIT);

    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = read_bits<unsigned char>(i * CHAR_BIT, CHAR_BIT);
    }

    return result;
}

gmprime::context_entropy gmprime::conditional_entropy(const size_t order, const size_t width) const {
    const auto values = symbols(width);

    if (values.size() <= order) {
        return {0, 0, 0, 0};
    }

    // Each symbol after the first order ones, with its context: counts of contexts and of (context, symbol)
    sequence_counts contexts(values, order), joint(values, order + 1);
    const size_t samples = values.size() - order;
    for (size_t start = 0; start < samples; ++start) {
        ++contexts[start];
        ++joint[start];
    }

    // H(X | C) = sum over (c, x) of p(c, x) log2(n(c) / n(c, x)), n(c) found under the same start
    double naive = 0;
    for (const auto &[start, count] : joint) {
        naive += static_cast<double>(count) * std::log2(static_cast<double>(contexts.at(start)) / static_cast<double>(count));
    }
    naive /= static_cast<double>(samples);

    // Miller-Madow, summed over the contexts: (symbols seen in the context - 1) / (2 n(c) ln 2) each,
    // weighted by n(c) / samples
    const auto correction = static_cast<double>(joint.size() - contexts.size()) / (2 * static_cast<double>(samples) * std::log(2.0));

    return {naive, naive + correction, contexts.size(), samples};
}

double gmprime::adaptive_code_length(const size_t order, const size_t width) const {
    const auto values = symbols(width);
    const double alphabet = std::ldexp(1.0, static_cast<int>(width));
    double bits = static_cast<double>(std::min(order, values.size()) * width);

    if (values.size() <= order) {
        return bits;
    }

    // Krichevsky-Trofimov: p(x | c) = (n(c, x) + 1/2) / (n(c) + alphabet / 2), counts updated after coding x
    sequence_counts contexts(values, order), joint(values, order + 1);
    for (size_t start = 0; start + order < values.size(); ++start) {
        auto &context_count = contexts[start];
        auto &joint_count = joint[start];
        bits -= std::log2((static_cast<double>(joint_count) + 0.5) / (static_cast<double>(context_count) + alphabet / 2));
        ++joint_count;
        ++context_count;
    }

    return bits;
}

std::vector<uint32_t> gmprime::longest_previous_factors() const {
    const auto text = bytes();
    const size_t n = text.size();

    if (n > std::numeric_limits<uint32_t>::max()) {
        throw std::length_error("Repetition analysis is limited to 4 GiB"); // LCOV_EXCL_LINE: needs a 4 GiB test
    }

    const auto sa = suffix_array(text);
    auto lcp = lcp_array(text, sa);
    std::vector<uint32_t> lpf(n, 0);

    // Crochemore, Ilie and Smyth: the longest previous factor of a suffix is its longest common prefix with
    // the nearest suffixes in sorted order that start before it, found with a stack in one pass. A virtual
    // suffix starting at -1 after the last one empties the stack.
    std::vector<size_t> stack;
    stack.push_back(0);
    for (size_t r = 1; r <= n; ++r) {
        uint32_t current = r < n ? lcp[r] : 0;

        while (!stack.empty()) {
            const auto top = stack.back();
            if (r == n || sa[r] < sa[top]) {
                lpf[sa[top]] = std::max(lcp[top], current);
                current = std::min(lcp[top], current);
            } else if (current <= lcp[top]) {
                lpf[sa[top]] = lcp[top];
            } else {
                break;
            }
            stack.pop_back();
        }

        if (r < n) {
            lcp[r] = current;
            stack.push_back(r);
        }
    }

    return lpf;
}

std::vector<size_t> gmprime::autocorrelation(const size_t max_distance) const {
    const size_t length = analyzed_bits();
    const gmprime magnitude = abs(*this);
    std::vector<size_t> result(max_distance + 1, 0);

    for (size_t d = 0; d <= max_distance && d < length; ++d) {
        const size_t overlap = length - d;
        // Bits differing from the one d further: the XOR of the number and its shift, without the top d
        // bits, where the shifted number is 0 and the XOR only repeats the number
        const size_t differing = (magnitude ^ (magnitude >> d)).popcount() - (magnitude >> overlap).popcount();
        result[d] = overlap - differing;
    }

    return result;
}

gmprime::repetition_stats gmprime::repetitions(const size_t min_length) const {
    const auto lpf = longest_previous_factors();
    const size_t threshold = std::max<size_t>(min_length, 1);
    repetition_stats stats{lpf.size(), 0, 0, 0};

    for (size_t i = 0; i < lpf.size();) {
        if (lpf[i] >= threshold) {
            ++stats.matches;
            stats.matched_bytes += lpf[i];
            i += lpf[i];
        } else {
            ++stats.literals;
            ++i;
        }
    }

    return stats;
}
