#ifndef GMPRIME_RECORDS_HPP
#define GMPRIME_RECORDS_HPP

// Tools for numeric data, alongside the coding toolkit: record widths and 16-bit images detected on bytes,
// integer linear prediction, and the lattice transform of floating point fields (decimal values stored in binary
// floating point, k * step for an integer k, replaced by k and a residual of a few units in the last place).
// Encoders analyze in floating point, but everything a decoder runs is integer arithmetic.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <optional>
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

} // GMPrime::coding

#endif //GMPRIME_RECORDS_HPP
