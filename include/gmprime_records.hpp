#ifndef GMPRIME_RECORDS_HPP
#define GMPRIME_RECORDS_HPP

// Tools for numeric data, alongside the coding toolkit: record widths and 16-bit images detected on bytes,
// integer linear prediction, and the lattice transform of floating point fields (decimal values stored in binary
// floating point, k * step for an integer k, replaced by k and a residual of a few units in the last place).
// Encoders analyze in floating point, but everything a decoder runs is integer arithmetic.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>
#include <gmprime_coding.hpp>

namespace GMPrime::coding {

// RECORDS

// Share of bytes equal to the byte d before, for d from 1 to max_distance (result[0] unused), over the first
// sample bytes, read every stride bytes: a fixed-size record of R bytes, or an image row, makes peaks at R and
// its multiples, at byte granularity where analyze_folds sees bits
inline std::vector<double> byte_autocorrelation(const std::vector<uint8_t>& data, const size_t max_distance, const size_t sample = 1 << 20,
                                                const size_t stride = 1) {
    const size_t n = std::min(data.size(), sample);
    std::vector<double> result(max_distance + 1, 0.0);
    for (size_t d = 1; d <= max_distance && d < n; ++d) {
        size_t equal = 0, count = 0;
        for (size_t i = d; i < n; i += stride) {
            equal += data[i] == data[i - d];
            ++count;
        }
        result[d] = count ? static_cast<double>(equal) / static_cast<double>(count) : 0.0;
    }
    return result;
}

struct record_analysis {
    size_t width = 0;    // record width in bytes, 0 if none stands out
    double score = 0;    // share of bytes equal to the one a record above
    double median = 0;   // median share over the distances searched
};

// The record width: the smallest distance whose byte equality is close to the best one (multiples of a record
// score about as well as the record), if it stands out from the median of the distances
inline record_analysis find_record_width(const std::vector<uint8_t>& data, const size_t max_width = 4096, const size_t sample = 1 << 18) {
    record_analysis result;
    const auto eq = byte_autocorrelation(data, max_width, sample);
    if (eq.size() < 3) {
        return result;
    }
    std::vector<double> sorted(eq.begin() + 1, eq.end());
    std::nth_element(sorted.begin(), sorted.begin() + static_cast<std::ptrdiff_t>(sorted.size() / 2), sorted.end());
    result.median = sorted[sorted.size() / 2];
    const double best = *std::max_element(eq.begin() + 1, eq.end());
    for (size_t d = 2; d < eq.size(); ++d) {
        // Distance 1 is a run, not a record
        if (eq[d] >= best - (best - result.median) * 0.1) {
            result.width = d;
            result.score = eq[d];
            break;
        }
    }
    return result;
}

// IMAGES OF 16-BIT PIXELS

struct image16_analysis {
    bool found = false;
    size_t start = 0;          // parity of the first byte of the pixels
    bool little_endian = true;
    size_t row_bytes = 0;      // 0 if no row width stands out
    double row_score = 1;      // mean difference with the sample a row above, relative to the samples around it
};

// Detects 16-bit samples (bytes alike two apart, not one apart; one byte of each pair, the high one, of lower
// entropy), their byte order (the one whose values change least from one sample to the next), and the row width
// of an image: the distance in samples where values differ least, searched from min_row to max_row samples
inline image16_analysis find_image16(const std::vector<uint8_t>& data, const size_t min_row = 16, const size_t max_row = 8192,
                                     const size_t sample = 1 << 19) {
    image16_analysis result;
    const size_t n = std::min(data.size(), sample) & ~size_t{1};
    if (n < 4 * max_row) {
        return result;
    }
    size_t equal1 = 0, equal2 = 0;
    std::vector<size_t> even(256, 0), odd(256, 0);
    for (size_t i = 2; i < n; ++i) {
        equal1 += data[i] == data[i - 1];
        equal2 += data[i] == data[i - 2];
        ++(i % 2 ? odd : even)[data[i]];
    }
    if (equal2 < n / 8 || equal2 < 2 * equal1) {
        return result;
    }
    const size_t high = gmprime::entropy(even) < gmprime::entropy(odd) ? 0 : 1;
    // Big-endian samples start at the high byte, little-endian ones just before it
    const auto pixels = [&](const size_t start, const bool little) {
        std::vector<int32_t> values((n - start) / 2);
        for (size_t j = 0; j < values.size(); ++j) {
            const size_t i = start + 2 * j;
            values[j] = little ? data[i] | data[i + 1] << 8 : data[i] << 8 | data[i + 1];
        }
        return values;
    };
    const auto roughness = [](const std::vector<int32_t>& v) {
        uint64_t sum = 0;
        for (size_t i = 1; i < v.size(); ++i) {
            sum += static_cast<uint64_t>(std::abs(v[i] - v[i - 1]));
        }
        return sum;
    };
    const auto big = pixels(high, false), little = pixels(1 - high, true);
    result.little_endian = roughness(little) < roughness(big);
    result.start = result.little_endian ? 1 - high : high;
    result.found = true;
    const auto& v = result.little_endian ? little : big;

    // Mean difference with the sample d before, each difference capped (a header or an edge would weigh on every
    // distance), then the row: the sharpest dip, the distance whose difference is lowest against those 4 samples
    // on either side (a smooth image differs least at small distances, its rows only stand out locally)
    std::vector<double> mad(max_row + 5, 0);
    const size_t stride = 3;
    const size_t top = std::min(max_row + 4, v.size() / 2);
    for (size_t d = min_row - 4; d <= top; ++d) {
        uint64_t sum = 0, count = 0;
        for (size_t i = d; i < v.size(); i += stride) {
            sum += static_cast<uint64_t>(std::min(std::abs(v[i] - v[i - d]), 1024));
            ++count;
        }
        mad[d] = static_cast<double>(sum) / static_cast<double>(count);
    }
    size_t best = 0;
    double best_ratio = 1;
    for (size_t d = min_row; d + 4 <= top; ++d) {
        const double around = (mad[d - 4] + mad[d + 4]) / 2;
        const double ratio = around > 0 ? mad[d] / around : 1;
        if (ratio < best_ratio) {
            best_ratio = ratio;
            best = d;
        }
    }
    result.row_score = best_ratio;
    if (best > 0 && best_ratio < 0.85) {
        result.row_bytes = 2 * best;
    }
    return result;
}

// NUMERIC PREDICTION

// Adaptive linear prediction of a value from features (neighbours as differences from a base value, each within
// +-2^24), with integer normalized LMS: weights in 16.16 fixed point, each step moving them along the error
// (capped at +-2^16) times the features over their energy, at rate 2^-rate_shift. Integers only, so that the
// decoder predicts exactly the same
class lms_predictor {
public:
    explicit lms_predictor(const size_t features, const int rate_shift = 4): _weights(features, 0), _x(features, 0), _shift(rate_shift) {
    }

    [[nodiscard]] size_t size() const {
        return _weights.size();
    }

    // The predicted difference from the base, for these features
    int64_t predict(const int64_t* features) {
        int64_t dot = 0;
        _energy = 64;
        for (size_t i = 0; i < _weights.size(); ++i) {
            _x[i] = features[i];
            dot += _weights[i] * features[i];
            _energy += features[i] * features[i];
        }
        _p = dot >> 16;
        return _p;
    }

    // Learns the actual difference
    void update(const int64_t actual) {
        const int64_t err = std::max<int64_t>(-65536, std::min<int64_t>(65536, actual - _p));
        for (size_t i = 0; i < _weights.size(); ++i) {
            _weights[i] += err * _x[i] * (int64_t{1} << (16 - _shift)) / _energy;
        }
    }

private:
    std::vector<int64_t> _weights, _x;
    int _shift;
    int64_t _energy = 64, _p = 0;
};

// A signed difference quantized on a logarithmic scale, about two levels per octave (0, then 1, 2, 3, 4, 5, 6...
// for 1, 2, 3, 4 to 5, 6 to 7, 8 to 11...): -127 to 127 over 64-bit differences. A prediction's distance from what
// is known of a value, as a context, shares its statistics between values far apart
inline int quantize_difference(const int64_t e) {
    const uint64_t a = e < 0 ? static_cast<uint64_t>(-(e + 1)) + 1 : static_cast<uint64_t>(e);
    if (a == 0) {
        return 0;
    }
    const int length = 64 - leading_zeros(a);
    const int level = length < 2 ? 1 : 2 * length - 2 + static_cast<int>((a >> (length - 2)) & 1);
    return e < 0 ? -level : level;
}

// Swaps the bytes of each pair from first on: 16-bit values of either byte order, high byte first
inline void swap_pairs(std::vector<uint8_t>& data, const size_t first) {
    for (size_t i = first; i + 1 < data.size(); i += 2) {
        std::swap(data[i], data[i + 1]);
    }
}

// The middle of the values left possible for a sample of bits bits once its first known bits are known (prefix,
// their value): the next bit is 1 when the sample is at least this middle
inline int64_t bit_middle(const int64_t prefix, const int known, const int bits) {
    return (prefix << (bits - known)) + (int64_t{1} << (bits - 1 - known));
}

// FLOAT LATTICES

// IEEE 754 formats: binary64 (double) or binary32 (float), in bytes of either order; or the IBM hexadecimal
// single precision of mainframes and seismic data (4 bytes, sign, exponent of 16 biased by 64, 24-bit fraction)
struct float_format {
    uint8_t size = 8;          // 8 or 4 bytes
    bool big_endian = false;
    bool ibm = false;          // IBM single precision (size 4); lattice steps are then powers of 2

    [[nodiscard]] int mantissa_bits() const { return size == 8 ? 52 : 23; }
    [[nodiscard]] int exponent_bits() const { return size == 8 ? 11 : 8; }
};

// The bits of a field read in its byte order, and written back
inline uint64_t read_field(const uint8_t* p, const float_format f) {
    uint64_t v = 0;
    for (int i = 0; i < f.size; ++i) {
        const int b = f.big_endian ? i : f.size - 1 - i;
        v = (v << 8) | p[b];
    }
    return v;
}

inline void write_field(uint8_t* p, const float_format f, uint64_t v) {
    for (int i = f.size - 1; i >= 0; --i) {
        const int b = f.big_endian ? i : f.size - 1 - i;
        p[b] = static_cast<uint8_t>(v & 0xFF);
        v >>= 8;
    }
}

// k * step rounded to the nearest representable value (ties to even), as the format's bit pattern, computed with
// integers only: exactly what IEEE multiplication gives on every platform, without depending on how a compiler
// evaluates floating point; +0 for k = 0. Nothing for a result that is not a finite normal number
// IBM single precision: k * 2^p exactly, normalized (first hexadecimal digit of the fraction not 0), or nothing
inline std::optional<uint64_t> ibm_lattice_bits(const int64_t k, const int p) {
    if (k == 0) {
        return uint64_t{0};
    }
    const uint64_t magnitude = k < 0 ? static_cast<uint64_t>(-(k + 1)) + 1 : static_cast<uint64_t>(k);
    const int length = 64 - leading_zeros(magnitude);
    // fraction = magnitude * 2^t, of 21 to 24 bits, with p + 24 - t a multiple of 4: value = fraction * 16^(E - 64) / 2^24
    int t = 24 - length;
    while (((p + 24 - t) % 4 + 4) % 4 != 0) {
        --t;
    }
    if (t < 0 && (magnitude & ((uint64_t{1} << -t) - 1)) != 0) {
        return std::nullopt;
    }
    const uint64_t fraction = t >= 0 ? magnitude << t : magnitude >> -t;
    const int exponent = 64 + (p + 24 - t) / 4;
    if (exponent < 0 || exponent > 127) {
        return std::nullopt;
    }
    return (k < 0 ? uint64_t{1} << 31 : 0) | static_cast<uint64_t>(exponent) << 24 | fraction;
}

// The power of 2 of an IBM step (a fraction of a single bit), or nothing
inline std::optional<int> ibm_step_power(const uint64_t step_bits) {
    const uint64_t fraction = step_bits & 0xFFFFFF;
    if (fraction == 0 || (fraction & (fraction - 1)) != 0 || (step_bits >> 31) != 0) {
        return std::nullopt;
    }
    return 4 * (static_cast<int>((step_bits >> 24) & 0x7F) - 64) - 24 + trailing_zeros(fraction);
}

// The 128-bit product of two 64-bit numbers, as its high and low words
struct wide_product {
    uint64_t high;
    uint64_t low;
};

inline wide_product multiply_wide(const uint64_t a, const uint64_t b) {
    const uint64_t a0 = a & 0xFFFFFFFF, a1 = a >> 32, b0 = b & 0xFFFFFFFF, b1 = b >> 32;
    const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    const uint64_t middle = (p00 >> 32) + (p01 & 0xFFFFFFFF) + (p10 & 0xFFFFFFFF);
    return {p11 + (p01 >> 32) + (p10 >> 32) + (middle >> 32), middle << 32 | (p00 & 0xFFFFFFFF)};
}

inline std::optional<uint64_t> lattice_bits(const int64_t k, const uint64_t step_bits, const float_format f) {
    if (f.ibm) {
        const auto p = ibm_step_power(step_bits);
        return p ? ibm_lattice_bits(k, *p) : std::nullopt;
    }
    const int mb = f.mantissa_bits(), eb = f.exponent_bits();
    const int bias = (1 << (eb - 1)) - 1;
    const uint64_t step_sign = step_bits >> (mb + eb);
    const int step_exp = static_cast<int>((step_bits >> mb) & ((uint64_t{1} << eb) - 1));
    uint64_t m = step_bits & ((uint64_t{1} << mb) - 1);
    if (step_exp == (1 << eb) - 1 || (step_exp == 0 && m == 0)) {
        return std::nullopt;
    }
    if (k == 0) {
        return uint64_t{0};
    }
    int e = step_exp == 0 ? 1 - bias - mb : step_exp - bias - mb;   // step = m * 2^e
    if (step_exp != 0) {
        m |= uint64_t{1} << mb;
    }
    const uint64_t magnitude = k < 0 ? static_cast<uint64_t>(-(k + 1)) + 1 : static_cast<uint64_t>(k);
    const auto p = multiply_wide(magnitude, m);
    const int length = p.high != 0 ? 128 - leading_zeros(p.high) : 64 - leading_zeros(p.low);
    // Bits [0, n) of the product, and whether any is set
    const auto low_bits_set = [&](const int n) {
        return n <= 0 ? false : n >= 64 ? p.low != 0 || (n > 64 && (p.high & (n >= 128 ? ~uint64_t{0} : (uint64_t{1} << (n - 64)) - 1)) != 0)
                                        : (p.low & ((uint64_t{1} << n) - 1)) != 0;
    };
    const auto bit = [&](const int i) {
        return static_cast<int>(i >= 64 ? (p.high >> (i - 64)) & 1 : (p.low >> i) & 1);
    };
    // Significand of mb + 1 bits: rounded to nearest, ties to even
    const int shift = length - (mb + 1);
    uint64_t significand;
    if (shift > 0) {
        significand = shift >= 64 ? p.high >> (shift - 64) : p.low >> shift | p.high << (64 - shift);
        if (bit(shift - 1) && (low_bits_set(shift - 1) || (significand & 1))) {
            ++significand;
            if (significand >> (mb + 1)) {
                significand >>= 1;
                ++e;
            }
        }
        e += shift;
    } else {
        significand = p.low << -shift;
        e += shift;
    }
    // value = significand * 2^e, significand in [2^mb, 2^(mb+1))
    const int biased = e + mb + bias;
    if (biased <= 0 || biased >= (1 << eb) - 1) {
        return std::nullopt;
    }
    const uint64_t sign = step_sign ^ (k < 0 ? 1 : 0);
    return sign << (mb + eb) | static_cast<uint64_t>(biased) << mb | (significand & ((uint64_t{1} << mb) - 1));
}

// A floating point field of each record, coded as an integer k on a lattice of step and a residual in units
// in the last place
struct lattice_field {
    size_t offset = 0;            // in the record
    float_format format;
    uint64_t step_bits = 0;       // the step, in the field's format
    int residual_bits = 4;        // residuals from -2^(residual_bits - 1) to 2^(residual_bits - 1) - 1 (0: none)
    bool delta = false;           // k coded as the difference with the previous record's
    bool offset_binary = false;   // k coded in offset binary (k + 2^(n - 1), in order: numeric contexts can predict
                                  // it), else zigzag (small magnitudes of either sign leave the high bytes 0)
};

// Values of a format as doubles, for the encoder's analysis only (the transform itself is exact on integers)
inline double field_value(const uint64_t bits, const float_format f) {
    if (f.ibm) {
        const double fraction = static_cast<double>(bits & 0xFFFFFF) / 16777216.0;
        const double value = std::ldexp(fraction, 4 * (static_cast<int>((bits >> 24) & 0x7F) - 64));
        return (bits >> 31) & 1 ? -value : value;
    }
    if (f.size == 8) {
        double d;
        std::memcpy(&d, &bits, 8);
        return d;
    }
    const auto b = static_cast<uint32_t>(bits);
    float x;
    std::memcpy(&x, &b, 4);
    return x;
}

inline uint64_t value_bits(const double v, const float_format f) {
    if (f.ibm) {
        // Steps only: the nearest power of 2
        const int p = static_cast<int>(std::lround(std::log2(std::fabs(v))));
        return ibm_lattice_bits(1, p).value_or(0);
    }
    if (f.size == 8) {
        uint64_t b;
        std::memcpy(&b, &v, 8);
        return b;
    }
    const auto x = static_cast<float>(v);
    uint32_t b;
    std::memcpy(&b, &x, 4);
    return b;
}

inline uint64_t zigzag(const int64_t v) {
    return v < 0 ? (static_cast<uint64_t>(-(v + 1)) << 1) | 1 : static_cast<uint64_t>(v) << 1;
}

inline int64_t unzigzag(const uint64_t z) {
    return (z & 1) ? -static_cast<int64_t>(z >> 1) - 1 : static_cast<int64_t>(z >> 1);
}

// The lattice representation of value bits v: k and the residual, if v is close enough to k * step
struct lattice_point {
    int64_t k = 0;
    int64_t residual = 0;
};

inline std::optional<lattice_point> to_lattice(const uint64_t v, const lattice_field& field, const double step) {
    const double x = field_value(v, field.format);
    if (!std::isfinite(x) || std::fabs(x / step) > 4e15) {
        return std::nullopt;
    }
    const auto k = static_cast<int64_t>(std::llround(x / step));
    const auto bits = lattice_bits(k, field.step_bits, field.format);
    if (!bits) {
        return std::nullopt;
    }
    const uint64_t sign_bit = uint64_t{1} << (field.format.size * 8 - 1);
    if ((v & sign_bit) != (*bits & sign_bit)) {
        return std::nullopt;
    }
    const int64_t residual = static_cast<int64_t>(v & (sign_bit - 1)) - static_cast<int64_t>(*bits & (sign_bit - 1));
    const int64_t limit = field.residual_bits > 0 ? int64_t{1} << (field.residual_bits - 1) : 0;
    if (field.residual_bits == 0 ? residual != 0 : residual < -limit || residual >= limit) {
        return std::nullopt;
    }
    return lattice_point{k, residual};
}

// The integer stored in the field's bytes for a coded k (k or its delta) and a residual, if it fits
inline std::optional<uint64_t> pack_point(const lattice_field& field, const int64_t coded, const int64_t residual) {
    const int total = field.format.size * 8;
    const int high = total - field.residual_bits;      // bits for the coded k
    const uint64_t bias = field.residual_bits > 0 ? uint64_t{1} << (field.residual_bits - 1) : 0;
    uint64_t c;
    if (field.offset_binary) {
        // Without residual bits, a 64-bit field holds any k
        const uint64_t half = uint64_t{1} << (high - 1);
        if (high < 64 && (coded < -static_cast<int64_t>(half) || coded >= static_cast<int64_t>(half))) {
            return std::nullopt;
        }
        c = static_cast<uint64_t>(coded) + half;
    } else {
        c = zigzag(coded);
        if ((c >> (high - 1)) != 0) {
            return std::nullopt;
        }
    }
    return c << field.residual_bits | static_cast<uint64_t>(residual + static_cast<int64_t>(bias));
}

inline lattice_point unpack_point(const lattice_field& field, const uint64_t u) {
    const int total = field.format.size * 8;
    const int high = total - field.residual_bits;
    const int64_t bias = field.residual_bits > 0 ? int64_t{1} << (field.residual_bits - 1) : 0;
    const uint64_t c = u >> field.residual_bits;
    const int64_t coded = field.offset_binary ? static_cast<int64_t>(c - (uint64_t{1} << (high - 1))) : unzigzag(c);
    return {coded, static_cast<int64_t>(u & ((uint64_t{1} << field.residual_bits) - 1)) - bias};
}

// Searches a lattice for a column of values (finite ones): the smallest gap between sorted distinct values,
// refined over growing spans so that k stays exact for large values. Encoder side, in floating point: the step
// found is stored, and only lattice_bits runs at decoding
inline std::optional<double> find_lattice_step(std::vector<double> values) {
    std::vector<double> magnitudes;
    for (const double v : values) {
        if (std::isfinite(v) && v != 0) {
            magnitudes.push_back(std::fabs(v));
        }
    }
    if (magnitudes.size() < 64) {
        return std::nullopt;
    }
    std::sort(magnitudes.begin(), magnitudes.end());
    magnitudes.erase(std::unique(magnitudes.begin(), magnitudes.end()), magnitudes.end());
    // The 99th percentile: stray huge values (a header read as a field) would hide every gap
    const double top = magnitudes[magnitudes.size() * 99 / 100];
    std::vector<double> gaps;
    for (size_t i = 1; i < magnitudes.size(); ++i) {
        const double g = magnitudes[i] - magnitudes[i - 1];
        if (g > top * 1e-12) {
            gaps.push_back(g);
        }
    }
    if (gaps.empty()) {
        return std::nullopt;
    }
    // A common small gap, the first percentile (the smallest one may join two stray values), or a fraction of it
    std::sort(gaps.begin(), gaps.end());
    // From the first percentile: a stray value (a header read as a field) would spoil every span
    const double low = magnitudes[magnitudes.size() / 100];
    for (int trial = 0; trial < 16; ++trial) {
        const int divisor = trial % 8 + 1;
        double step = (trial < 8 ? gaps.front() : gaps[gaps.size() / 100]) / divisor;
        for (const double v : magnitudes) {
            if (v < low || v > top) {
                continue;
            }
            const double steps = (v - low) / step;
            const double m = std::round(steps);
            if (m > 0 && std::fabs(steps - m) < 0.05) {
                step = (v - low) / m;
            }
        }
        // The values must be multiples of the step, not only their differences
        size_t near = 0;
        for (const double v : magnitudes) {
            const double q = v / step;
            near += std::fabs(q - std::round(q)) < 0.01;
        }
        if (near >= magnitudes.size() * 9 / 10) {
            return step;
        }
    }
    return std::nullopt;
}

// Bits the transformed field would take, estimated over a sample: each of its bytes coded by an adaptive model,
// the better of order 0 and of the same byte in the record above as context (what column contexts reach at
// least: high bytes that never change cost nothing, slowly changing ones little); an escaped field costs its size
inline double lattice_cost(const std::vector<uint8_t>& data, const size_t width, const lattice_field& field, const size_t records, size_t* escapes) {
    const double step = field_value(field.step_bits, field.format);
    const size_t size = field.format.size;
    // Counts with an extra total per context: [byte position][context][symbol or 256]
    std::vector<uint32_t> order0(size * 257, 0), above(size * 256 * 257, 0);
    std::vector<double> cost0(size, 0), cost1(size, 0);
    int64_t previous = 0;
    uint64_t last = 0;
    *escapes = 0;
    const int total = field.format.size * 8;
    const double alpha = 1.0 / 16;
    for (size_t r = 0; r < records; ++r) {
        const auto v = read_field(&data[r * width + field.offset], field.format);
        const auto point = to_lattice(v, field, step);
        const int64_t coded = point ? (field.delta ? point->k - previous : point->k) : 0;
        const auto packed = point ? pack_point(field, coded, point->residual) : std::nullopt;
        if (!packed) {
            ++*escapes;
            continue;
        }
        previous = point->k;
        const uint64_t u = *packed;
        for (size_t i = 0; i < size; ++i) {
            const size_t b = (u >> (8 * i)) & 0xFF, a = (last >> (8 * i)) & 0xFF;
            uint32_t* c0 = &order0[i * 257];
            uint32_t* c1 = &above[(i * 256 + a) * 257];
            cost0[i] -= std::log2((c0[b] + alpha) / (c0[256] + 256 * alpha));
            cost1[i] -= std::log2((c1[b] + alpha) / (c1[256] + 256 * alpha));
            ++c0[b], ++c0[256], ++c1[b], ++c1[256];
        }
        last = u;
    }
    double bits = static_cast<double>(*escapes) * total;
    for (size_t i = 0; i < size; ++i) {
        bits += std::min(cost0[i], cost1[i]);
    }
    return bits;
}

// Finds the floating point fields of records of width bytes that sit on a lattice, best savings first, without
// overlaps. Encoder side
inline std::vector<lattice_field> find_lattice_fields(const std::vector<uint8_t>& data, const size_t width, const size_t max_records = 20000) {
    std::vector<lattice_field> found;
    if (width < 4) {
        return found;
    }
    const size_t records = std::min(data.size() / width, max_records);
    if (records < 256) {
        return found;
    }
    struct candidate {
        lattice_field field;
        double saving;
    };
    std::vector<candidate> candidates;
    const float_format formats[] = {{8, false, false}, {8, true, false}, {4, false, false}, {4, true, false}, {4, true, true}};
    for (const auto f : formats) {
        const uint8_t size = f.size;
        for (size_t offset = 0; offset + size <= width; ++offset) {
            std::vector<double> values;
            for (size_t r = 0; r < records; ++r) {
                values.push_back(field_value(read_field(&data[r * width + offset], f), f));
            }
            const auto step = find_lattice_step(values);
            if (!step) {
                continue;
            }
            lattice_field field;
            field.offset = offset;
            field.format = f;
            field.step_bits = value_bits(*step, f);
            // A step at the precision of the values themselves (any float in [1, 2) is a multiple of 2^-23)
            // only renames the mantissa: the lattice must hold the values in clearly fewer bits
            std::vector<double> sorted_k;
            for (const double v : values) {
                if (std::isfinite(v)) {
                    sorted_k.push_back(std::fabs(v) / *step);
                }
            }
            std::nth_element(sorted_k.begin(), sorted_k.begin() + static_cast<std::ptrdiff_t>(sorted_k.size() * 99 / 100), sorted_k.end());
            if (sorted_k.empty() || sorted_k[sorted_k.size() * 99 / 100] > std::ldexp(1.0, f.mantissa_bits() - 2)) {
                continue;
            }
            double best = 0;
            lattice_field chosen;
            const std::vector<int> residuals = f.ibm ? std::vector<int>{0} : std::vector<int>{0, 2, 3, 4, 5};
            for (const int residual_bits : residuals) {
                for (const int variant : {0, 1, 2, 3}) {
                    field.residual_bits = residual_bits;
                    field.delta = variant & 1;
                    field.offset_binary = variant & 2;
                    size_t escapes = 0;
                    const double cost = lattice_cost(data, width, field, records, &escapes);
                    // Escaped fields cost a flag and their bytes; the raw field costs its size, a rough bound
                    const double saving = static_cast<double>(records) * size * 8 - cost;
                    if (escapes < records / 4 && saving > best) {
                        best = saving;
                        chosen = field;
                    }
                }
            }
            if (best > static_cast<double>(records) * 8) {
                candidates.push_back({chosen, best});
            }
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const candidate& a, const candidate& b) { return a.saving > b.saving; });
    std::vector<bool> used(width, false);
    for (const auto& c : candidates) {
        if (found.size() == 8) {
            break;
        }
        bool overlap = false;
        for (size_t i = 0; i < c.field.format.size; ++i) {
            overlap = overlap || used[c.field.offset + i];
        }
        if (!overlap) {
            for (size_t i = 0; i < c.field.format.size; ++i) {
                used[c.field.offset + i] = true;
            }
            found.push_back(c.field);
        }
    }
    std::sort(found.begin(), found.end(), [](const lattice_field& a, const lattice_field& b) { return a.offset < b.offset; });
    return found;
}

// The transform: each record of width bytes gets a first byte of flags (bit j: field j kept as it was), and each
// field on its lattice is replaced by zigzag(k or its delta) << residual_bits | residual biased, big-endian in
// the field's bytes: the high bytes, mostly zero or alike, come first, and the field's information sits in its
// low bytes. Bytes after the last whole record are kept. Records grow by one byte
inline std::vector<uint8_t> lattice_encode(const std::vector<uint8_t>& data, const size_t width, const std::vector<lattice_field>& fields) {
    std::vector<uint8_t> out;
    const size_t records = width ? data.size() / width : 0;
    out.reserve(data.size() + records + 16);
    std::vector<int64_t> previous(fields.size(), 0);
    std::vector<double> steps;
    for (const auto& field : fields) {
        steps.push_back(field_value(field.step_bits, field.format));
    }
    for (size_t r = 0; r < records; ++r) {
        const uint8_t* rec = &data[r * width];
        const size_t start = out.size();
        out.push_back(0);
        out.insert(out.end(), rec, rec + width);
        uint8_t flags = 0;
        for (size_t j = 0; j < fields.size(); ++j) {
            const auto& field = fields[j];
            const auto v = read_field(rec + field.offset, field.format);
            const auto point = to_lattice(v, field, steps[j]);
            const int64_t coded = point ? (field.delta ? point->k - previous[j] : point->k) : 0;
            const auto packed = point ? pack_point(field, coded, point->residual) : std::nullopt;
            if (!packed) {
                flags |= static_cast<uint8_t>(1 << j);
                continue;
            }
            previous[j] = point->k;
            const uint64_t u = *packed;
            for (int i = 0; i < field.format.size; ++i) {
                out[start + 1 + field.offset + static_cast<size_t>(i)] = static_cast<uint8_t>(u >> (8 * (field.format.size - 1 - i)));
            }
        }
        out[start] = flags;
    }
    out.insert(out.end(), data.begin() + static_cast<std::ptrdiff_t>(records * width), data.end());
    return out;
}

// The inverse: nothing if the encoded data is not a transform of records of this width and fields
inline std::optional<std::vector<uint8_t>> lattice_decode(const std::vector<uint8_t>& encoded, const size_t width, const std::vector<lattice_field>& fields) {
    std::vector<uint8_t> out;
    const size_t records = encoded.size() / (width + 1);
    out.reserve(encoded.size());
    std::vector<int64_t> previous(fields.size(), 0);
    for (size_t r = 0; r < records; ++r) {
        const uint8_t* rec = &encoded[r * (width + 1)];
        const uint8_t flags = rec[0];
        const size_t start = out.size();
        out.insert(out.end(), rec + 1, rec + 1 + width);
        for (size_t j = 0; j < fields.size(); ++j) {
            if (flags & (1 << j)) {
                continue;
            }
            const auto& field = fields[j];
            uint64_t u = 0;
            for (int i = 0; i < field.format.size; ++i) {
                u = (u << 8) | rec[1 + field.offset + static_cast<size_t>(i)];
            }
            const auto point = unpack_point(field, u);
            const int64_t residual = point.residual, coded = point.k;
            const int64_t k = field.delta ? previous[j] + coded : coded;
            previous[j] = k;
            const auto bits = lattice_bits(k, field.step_bits, field.format);
            if (!bits) {
                return std::nullopt;
            }
            const uint64_t sign_bit = uint64_t{1} << (field.format.size * 8 - 1);
            const uint64_t v = (*bits & sign_bit) | static_cast<uint64_t>(static_cast<int64_t>(*bits & (sign_bit - 1)) + residual);
            write_field(&out[start + field.offset], field.format, v);
        }
    }
    // Records are width + 1 bytes: the tail is what remains after them
    const size_t tail = records * (width + 1);
    out.insert(out.end(), encoded.begin() + static_cast<std::ptrdiff_t>(tail), encoded.end());
    return out;
}

// NUMBERS IN TEXT. A reversible transform replaces each decimal number of a text by a token: a marker byte (absent
// from the text), its layout, its digits in a balanced 4-bit code; numeric_text_reader reads that stream causally,
// for the models that code it. A model then predicts each number as a number: rescale_decimal brings predictions
// to its precision, decimal_position tells where they stand against the digits still possible (the decimal
// counterpart of bit_middle), split_at what each coded bit of a digit splits. Programs print binary floats with
// more digits than they hold, and dates computed from counts of seconds: float32_print and civil_from_seconds
// compute those, with integers only.
//
// A number is a run of 1 to 18 digits, with a '-' sign when the '-' does not follow a letter, a digit or a '.',
// and a fraction ('.' and digits, 18 digits in all) unless it belongs to a dotted sequence (an IP address, a
// version, a date like 2005.06.03): its layout is the sign, the count of integer digits L (leading zeros
// included) and of fraction digits F, its value the integer m of its L + F digits. The layout byte holds
// sign << 7 | (F * 18 + L - 1) for F up to 6, else sign << 7 | 127 followed by L and F; then come the L + F
// digits, without the sign and the point. Every quirk of writing (leading zeros, trailing zeros) is in the layout,
// so that the text comes back exactly. Knowing the layout first, a model knows the weight of every digit. Anything
// else stays text between numbers ('+', the 'e-' of 1e-5, the ',' of 1,234, the '.' of .5); a run of more than
// 18 digits is cut in numbers of 18.

struct text_number {
    bool negative = false;
    int int_digits = 1;      // L, leading zeros included, 1 to 18
    int frac_digits = 0;     // F, 0 to 17
    uint64_t mantissa = 0;   // the L + F digits as one integer
};

inline uint64_t power10(const int n) {
    static const auto table = [] {
        std::array<uint64_t, 20> t{};
        t[0] = 1;
        for (size_t i = 1; i < t.size(); ++i) {
            t[i] = t[i - 1] * 10;
        }
        return t;
    }();
    return table[static_cast<size_t>(n)];
}

// A byte absent from the text, for the marker: control bytes first
inline std::optional<uint8_t> numeric_text_marker(const std::vector<uint8_t>& text) {
    std::array<bool, 256> seen{};
    for (const auto b : text) {
        seen[b] = true;
    }
    static const uint8_t preferred[] = {0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F, 0x03, 0x04, 0x05, 0x06, 0x07, 0x0E, 0x0F, 0x10, 0x11, 0x12, 0x13,
                                        0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x7F, 0x08, 0x0B, 0x0C};
    for (const auto b : preferred) {
        if (!seen[b]) {
            return b;
        }
    }
    for (int b = 0x80; b < 0x100; ++b) {
        if (!seen[static_cast<size_t>(b)]) {
            return static_cast<uint8_t>(b);
        }
    }
    return std::nullopt;
}

namespace detail {

inline bool is_digit(const uint8_t c) {
    return c >= '0' && c <= '9';
}

inline bool is_alnum(const uint8_t c) {
    return is_digit(c) || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

} // namespace detail

// Digits in the stream: 0x30 | a 4-bit code, monotone, whose binary tree splits the digits evenly (0-4 | 5-9, then
// 0-1 | 2-4 and 5-6 | 7-9, ...): each coded bit halves the values left, as the bits of a binary number do
inline uint8_t digit_code(const int d) {
    static const uint8_t codes[10] = {0, 2, 4, 6, 7, 8, 10, 12, 14, 15};
    return static_cast<uint8_t>(0x30 | codes[d]);
}

// The digit of a code, or -1
inline int code_digit(const uint8_t b) {
    static const int8_t digits[16] = {0, -1, 1, -1, 2, -1, 3, 4, 5, -1, 6, -1, 7, -1, 8, 9};
    return (b & 0xF0) == 0x30 ? digits[b & 15] : -1;
}

// At a node of a digit's code (its known low bits after a leading 1, 1 to 15): the digits left, from low to high,
// and the first one of them whose next bit is 1 (split > high when none is)
struct digit_split {
    int low = 0, split = 0, high = 9;
};

inline digit_split split_at(const uint32_t node) {
    static const auto table = [] {
        std::array<digit_split, 16> t{};
        for (uint32_t n = 1; n < 16; ++n) {
            const int known = n >= 8 ? 3 : n >= 4 ? 2 : n >= 2 ? 1 : 0;
            const uint32_t prefix = n & ((1u << known) - 1);
            digit_split s{10, 10, -1};
            for (int d = 0; d < 10; ++d) {
                const uint32_t c = digit_code(d) & 15u;
                if ((c >> (4 - known)) != prefix) {
                    continue;
                }
                s.low = std::min(s.low, d);
                s.high = std::max(s.high, d);
                if (((c >> (3 - known)) & 1) && d < s.split) {
                    s.split = d;
                }
            }
            t[n] = s;
        }
        return t;
    }();
    return table[node & 15];
}

// Layout bytes of a number
inline void put_layout(std::vector<uint8_t>& out, const text_number& n) {
    const int sign = n.negative ? 0x80 : 0;
    if (n.frac_digits <= 6) {
        out.push_back(static_cast<uint8_t>(sign | (n.frac_digits * 18 + n.int_digits - 1)));
    } else {
        out.push_back(static_cast<uint8_t>(sign | 127));
        out.push_back(static_cast<uint8_t>(n.int_digits));
        out.push_back(static_cast<uint8_t>(n.frac_digits));
    }
}

inline std::vector<uint8_t> numeric_text_encode(const std::vector<uint8_t>& text, const uint8_t marker) {
    std::vector<uint8_t> out;
    const size_t n = text.size();
    // Allocation failures abort: no cleanup to plan for
    [&]() noexcept {
        out.reserve(n + n / 4);
        size_t i = 0;
        while (i < n) {
            const uint8_t c = text[i];
            const bool sign = c == '-' && i + 1 < n && detail::is_digit(text[i + 1]) && (i == 0 || (!detail::is_alnum(text[i - 1]) && text[i - 1] != '.'));
            if (!detail::is_digit(c) && !sign) {
                out.push_back(c);
                ++i;
                continue;
            }
            text_number num;
            size_t p = sign ? i + 1 : i;
            num.negative = sign;
            const size_t start = p;
            uint64_t m = 0;
            int digits = 0;
            while (p < n && detail::is_digit(text[p]) && digits < 18) {
                m = m * 10 + (text[p] - '0');
                ++digits;
                ++p;
            }
            num.int_digits = digits;
            // A fraction, unless the number is part of a dotted sequence
            const bool dotted_before = start >= 2 && text[start - 1] == '.' && detail::is_digit(text[start - 2]);
            if (!dotted_before && digits < 18 && p + 1 < n && text[p] == '.' && detail::is_digit(text[p + 1])) {
                size_t q = p + 1;
                while (q < n && detail::is_digit(text[q])) {
                    ++q;
                }
                const int frac = static_cast<int>(q - p - 1);
                const bool dotted_after = q + 1 < n && text[q] == '.' && detail::is_digit(text[q + 1]);
                if (!dotted_after && digits + frac <= 18) {
                    for (size_t k = p + 1; k < q; ++k) {
                        m = m * 10 + (text[k] - '0');
                    }
                    num.frac_digits = frac;
                    p = q;
                }
            }
            num.mantissa = m;
            out.push_back(marker);
            put_layout(out, num);
            for (size_t k = start; k < p; ++k) {
                if (text[k] != '.') {
                    out.push_back(digit_code(text[k] - '0'));
                }
            }
            i = p;
        }
    }();
    return out;
}

// The text, or nothing if the stream is not a valid transform
inline std::optional<std::vector<uint8_t>> numeric_text_decode(const std::vector<uint8_t>& encoded, const uint8_t marker) {
    std::vector<uint8_t> out;
    out.reserve(encoded.size() + encoded.size() / 2);
    const size_t n = encoded.size();
    size_t i = 0;
    while (i < n) {
        if (encoded[i] != marker) {
            out.push_back(encoded[i++]);
            continue;
        }
        if (++i >= n) {
            return std::nullopt;
        }
        text_number num;
        const uint8_t layout = encoded[i++];
        num.negative = (layout & 0x80) != 0;
        const int code = layout & 0x7F;
        if (code == 127) {
            if (i + 2 > n) {
                return std::nullopt;
            }
            num.int_digits = encoded[i++];
            num.frac_digits = encoded[i++];
        } else if (code < 126) {
            num.int_digits = code % 18 + 1;
            num.frac_digits = code / 18;
        } else {
            return std::nullopt;
        }
        const int digits = num.int_digits + num.frac_digits;
        if (num.int_digits < 1 || digits > 18) {
            return std::nullopt;
        }
        if (i + static_cast<size_t>(digits) > n) {
            return std::nullopt;
        }
        if (num.negative) {
            out.push_back('-');
        }
        for (int d = 0; d < digits; ++d) {
            const int digit = code_digit(encoded[i++]);
            if (digit < 0) {
                return std::nullopt;
            }
            if (d == num.int_digits) {
                out.push_back('.');
            }
            out.push_back(static_cast<uint8_t>('0' + digit));
        }
    }
    return out;
}

// Share of digits in the text, in 1/10000
inline size_t digit_density(const std::vector<uint8_t>& text, const size_t sample = 65536) {
    const size_t n = std::min(text.size(), sample);
    size_t digits = 0;
    for (size_t i = 0; i < n; ++i) {
        digits += detail::is_digit(text[i]);
    }
    return n ? digits * 10000 / n : 0;
}

// DECIMAL PREDICTION. A value v with f fraction digits brought to t fraction digits (rounded half away from zero
// when digits are dropped), within +-10^18
inline int64_t rescale_decimal(const int64_t v, const int from, const int to) {
    if (to >= from) {
        const auto p = static_cast<int64_t>(power10(std::min(to - from, 18)));
        const int64_t limit = static_cast<int64_t>(power10(18)) / p;
        return std::max(-limit, std::min(limit, v)) * p;
    }
    const auto p = static_cast<int64_t>(power10(std::min(from - to, 18)));
    return v >= 0 ? (v + p / 2) / p : -((-v + p / 2) / p);
}

// The decimal counterpart of bit_middle and quantize_difference, as a context: where a prediction of a value of
// digits stands once its first digits are known (low: their value followed by zeros, unit: the weight of the next
// digit). Below or above the values left (64 + or 80 + the distance in units, on a log scale), or the digit it
// falls in, times 4, plus its quarter of that digit, so that carries are seen coming (0 to 39)
inline uint64_t decimal_position(const int64_t prediction, const int64_t low, const int64_t unit) {
    const int64_t rel = prediction - low;
    const auto lg = [](const int64_t d) { return static_cast<uint64_t>(std::min(15, 64 - leading_zeros(static_cast<uint64_t>(d)))); };
    if (rel < 0) {
        return 64 + lg((-rel + unit - 1) / unit);
    }
    if (rel >= 10 * unit) {
        return 80 + lg((rel - 10 * unit) / unit + 1);
    }
    const int64_t digit = rel / unit;
    const int64_t quarter = unit >= 4 ? (rel % unit) * 4 / unit : 0;
    return static_cast<uint64_t>(digit * 4 + quarter);
}

// BINARY FLOATS WRITTEN IN DECIMAL. Programs print binary floats: a float32 printed with more digits than it
// holds (127.49f as 127.489998, 35.99583f as 35.995834350586) leaves digits that its value computes. Integers
// only (128-bit products), so that a decoder predicts exactly the same

// m * 2^e, m below 2^24 (a float32 magnitude)
struct binary_float {
    uint64_t m = 0;
    int e = 0;
};

namespace detail {

inline int bit_length(const uint64_t v) {
    return v == 0 ? 0 : 64 - leading_zeros(v);
}

// (high, low) / d for high < d, so that the quotient fits 64 bits: the quotient and the remainder (long division
// by 32-bit digits, as in Hacker's Delight)
inline std::pair<uint64_t, uint64_t> divide_wide(const uint64_t high, const uint64_t low, uint64_t d) {
    constexpr uint64_t base = uint64_t{1} << 32;
    const int s = leading_zeros(d);
    d <<= s;
    const uint64_t d1 = d >> 32, d0 = d & 0xFFFFFFFF;
    const uint64_t n32 = s == 0 ? high : high << s | low >> (64 - s);
    const uint64_t n10 = low << s, n1 = n10 >> 32, n0 = n10 & 0xFFFFFFFF;
    uint64_t q1 = n32 / d1, r = n32 - q1 * d1;
    while (q1 >= base || q1 * d0 > (r << 32 | n1)) {
        --q1;
        r += d1;
        if (r >= base) {
            break;
        }
    }
    const uint64_t n21 = (n32 << 32 | n1) - q1 * d;
    uint64_t q0 = n21 / d1;
    r = n21 - q0 * d1;
    while (q0 >= base || q0 * d0 > (r << 32 | n0)) {
        --q0;
        r += d1;
        if (r >= base) {
            break;
        }
    }
    return {q1 << 32 | q0, ((n21 << 32 | n0) - q0 * d) >> s};
}

// The quotient q of (high, low) / d rounded half to even, given the remainder
inline uint64_t round_quotient(const uint64_t q, const uint64_t remainder, const uint64_t d) {
    const uint64_t half = d - remainder;   // remainder * 2 against d, without overflow
    return remainder > half || (remainder == half && (q & 1)) ? q + 1 : q;
}

} // namespace detail

// The float32 nearest to n / 10^f, f from 0 to 18: always a normal float32, n / 10^f being between 10^-18 and
// 2^64. The exponent is estimated from the bit lengths (within 2), then corrected
inline binary_float nearest_float32(const uint64_t n, const int f) {
    if (n == 0) {
        return {};
    }
    const int nb = detail::bit_length(n);
    const uint64_t p = power10(f);
    // 3402 / 1024 is above log2(10): the estimate is never too low, the quotient starts below 2^24
    int e = nb - 24 - ((f * 3402) >> 10);
    while (true) {
        // m = round(n * 2^-e / 10^f) in [2^23, 2^24); the quotient stays below 2^26, within 64 bits
        uint64_t q, remainder, d = p;
        if (e <= 0) {
            const int shift = -e;
            const uint64_t high = shift == 0 ? 0 : shift >= 64 ? n << (shift - 64) : n >> (64 - shift);
            const uint64_t low = shift >= 64 ? 0 : n << shift;
            std::tie(q, remainder) = detail::divide_wide(high, low, p);
        } else {
            d = p << e;
            q = n / d;
            remainder = n % d;
        }
        if (q < (1u << 23)) {
            --e;
        } else {
            binary_float r{detail::round_quotient(q, remainder, d), e};
            if (r.m == (1u << 24)) {
                r.m >>= 1;
                ++r.e;
            }
            return r;
        }
    }
}

// m * 2^e written with f decimals, as an integer of digits (round half to even), if below 10^19
inline std::optional<uint64_t> render_decimal(const binary_float x, const int f) {
    const auto scaled = multiply_wide(x.m, power10(f));
    if (x.e >= 0) {
        // Shifted left, it must stay below 10^19, so within 64 bits
        if (scaled.high != 0 || detail::bit_length(scaled.low) + x.e > 64) {
            return std::nullopt;
        }
        const uint64_t v = scaled.low << x.e;
        return v < power10(19) ? std::optional<uint64_t>(v) : std::nullopt;
    }
    // Shifted right by s, rounded half to even: the bits shifted out decide, against half of 2^s
    const int s = -x.e;
    if (s >= 128) {
        return uint64_t{0};
    }
    const auto bit = [&](const int i) { return i >= 64 ? (scaled.high >> (i - 64)) & 1 : (scaled.low >> i) & 1; };
    const auto below = [&](const int n) {   // any of bits [0, n) set
        return n <= 0 ? false : n >= 64 ? scaled.low != 0 || (n > 64 && (scaled.high & ((n >= 128 ? 0 : uint64_t{1} << (n - 64)) - 1)) != 0)
                                        : (scaled.low & ((uint64_t{1} << n) - 1)) != 0;
    };
    if (s < 64 && scaled.high >> s != 0) {
        return std::nullopt;    // 64 bits or more
    }
    uint64_t v = s >= 64 ? scaled.high >> (s - 64) : scaled.low >> s | (s == 0 ? 0 : scaled.high << (64 - s));
    if (bit(s - 1) && (below(s - 1) || (v & 1))) {
        ++v;
    }
    return v < power10(19) ? std::optional<uint64_t>(v) : std::nullopt;
}

// What a float32 near n / 10^f prints with f decimals (f from 0 to 18): n itself if n is such a print
inline std::optional<uint64_t> float32_print(const uint64_t n, const int f) {
    return render_decimal(nearest_float32(n, f), f);
}

// TIMES. A count of seconds since 1970 written in a text is often written again as a date and a time: the civil
// date of a day count (proleptic Gregorian calendar, integers only)
struct civil_time {
    int64_t year = 1970;
    int month = 1, day = 1, hour = 0, minute = 0, second = 0;
};

inline civil_time civil_from_seconds(const int64_t seconds) {
    const int64_t days = seconds >= 0 ? seconds / 86400 : -((-seconds + 86399) / 86400);
    const int64_t rest = seconds - days * 86400;
    const int64_t z = days + 719468;
    const int64_t era = (z >= 0 ? z : z - 146096) / 146097;
    const int64_t doe = z - era * 146097;
    const int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const int64_t mp = (5 * doy + 2) / 153;
    civil_time t;
    t.day = static_cast<int>(doy - (153 * mp + 2) / 5 + 1);
    t.month = static_cast<int>(mp < 10 ? mp + 3 : mp - 9);
    t.year = yoe + era * 400 + (t.month <= 2);
    t.hour = static_cast<int>(rest / 3600);
    t.minute = static_cast<int>(rest / 60 % 60);
    t.second = static_cast<int>(rest % 60);
    return t;
}

// The causal reader: fed the transformed stream byte by byte, it says where the next byte stands (text, layout,
// or a byte of a number's value) and gives each number once its last byte is known
class numeric_text_reader {
public:
    enum state_kind { text, layout, layout_l, layout_f, value };

    explicit numeric_text_reader(const uint8_t marker): _marker(marker) {
    }

    [[nodiscard]] state_kind state() const {
        return _state;
    }
    [[nodiscard]] const text_number& number() const {
        return _number;
    }
    // While in a value: its digits, the index of the next one, the value of the digits before it
    [[nodiscard]] int value_digits() const {
        return _bytes;
    }
    [[nodiscard]] int value_index() const {
        return _index;
    }
    [[nodiscard]] uint64_t value_prefix() const {
        return _prefix;
    }

    // Returns true when the byte completed a number
    bool byte_added(const uint8_t byte) {
        switch (_state) {
        case text:
            if (byte == _marker) {
                _state = layout;
                _number = text_number{};
            }
            return false;
        case layout: {
            _number.negative = (byte & 0x80) != 0;
            const int code = byte & 0x7F;
            if (code == 127) {
                _state = layout_l;
                return false;
            }
            _number.int_digits = code % 18 + 1;
            _number.frac_digits = std::min(code / 18, 17);
            start_value();
            return false;
        }
        case layout_l:
            _number.int_digits = std::max(1, std::min<int>(byte, 18));
            _state = layout_f;
            return false;
        case layout_f:
            _number.frac_digits = std::max(0, std::min<int>(byte, 18 - _number.int_digits));
            start_value();
            return false;
        default:    // value
            _prefix = _prefix * 10 + static_cast<uint64_t>(std::max(0, code_digit(byte)));
            if (++_index == _bytes) {
                _number.mantissa = _prefix;
                _state = text;
                return true;
            }
            return false;
        }
    }

private:
    void start_value() {
        _bytes = _number.int_digits + _number.frac_digits;
        _index = 0;
        _prefix = 0;
        _state = value;
    }

    uint8_t _marker;
    state_kind _state = text;
    text_number _number;
    int _bytes = 0, _index = 0;
    uint64_t _prefix = 0;
};

} // GMPrime::coding

#endif //GMPRIME_RECORDS_HPP
