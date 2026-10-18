// gmz: a context-mixing compressor built on gmprime and its coding toolkit.
//
//   gmz [options] c INPUT OUTPUT    compress
//   gmz [options] d INPUT OUTPUT    decompress
//   gmz [options] t INPUT           compress, decompress, check and report, without writing anything
//
//   -m MODE   auto (default: the structure plans compared), text, fold (full and bitmap variants tried), fold-full,
//             bitmap
//   -C MODE   capitals transform for text: auto (tried both ways, default), on, off
//   -S BYTES  ways compared on the first BYTES bytes, then the data compressed once (default 65536); 0 compresses
//             the data both ways and keeps the smaller message, slower and sometimes a little smaller
//   -T BITS   context table size in bits of slot count, instead of one sized for the data (the decompressor
//             needs the same value)
//   -1, -2    faster levels: text compresses 2 to 2.5 times faster and decompresses 1.6 to 1.9 times faster, for
//             messages 1.1% and 1.8% larger over the Canterbury corpus. Fewer contexts and mixers, nibble buckets,
//             smaller tables, no capitals comparison; folded data only loses its text contexts; -2 also drops the
//             chained maps. The level is stored in the message for the decompressor
//   -x MODE   machine code: auto (default), off, bin (sparse contexts for binary data), 32 or 64 (x86 code)
//
// The file is first analyzed with gmprime::analyze_folds. Data that folds (bitmap rows, fixed-size records)
// is coded with bit and byte templates across the fold, in 2D and 3D; anything else with the text model:
// context orders, words, line structure, last occurrences, indirect modeling, four mixers and chained maps.
// Binary data (10% of control bytes) adds sparse contexts; x86 code (frequent CALL and JMP with near targets)
// has its branch targets made absolute, and instruction contexts from a light decoder in place of the text
// contexts.
// In auto mode, other data is analyzed for structure, and each way it suggests compresses the first bytes, the
// smallest message deciding: text, the fold of analyze_folds whatever its strength, fixed-size records (the byte
// distance where bytes repeat most, with column contexts), their floating point fields on a decimal lattice
// (k * step stored as k and a residual of a few units in the last place; a 4-byte field also as a signal,
// predicted from the records before), and images of 16-bit samples (byte order and row width found, pixels
// coded high byte first, with numeric contexts: where predictions stand from the values the known bits allow).
// The compressed message is a gmprime: the arithmetic coder's output, saved and loaded like any number.
#include <chrono>
#include <cstdio>
#include <memory>
#include <optional>
#include <string>
#include <gmprime_coding.hpp>
#include <gmprime_records.hpp>

using namespace GMPrime;
using namespace GMPrime::coding;

namespace {

// What the model uses, decided by the analysis and stored in the message header
struct config {
    bool fold = false;      // folded model: templates across the fold
    bool bitmap = false;    // folded model without byte contexts, for pixels
    size_t width = 0;       // fold width in bits
    size_t layer = 0;       // rows per layer of the 3D view, 0 for none
    size_t record = 0;      // record width in bytes for the column contexts, 0 for none
    size_t lattice_width = 0;                // records of the lattice transform, 0 for none
    std::vector<lattice_field> lattice;      // its fields
    size_t sample_bytes = 0;  // numeric samples of 2 or 4 bytes, high byte first, 0 for none
    size_t sample_stride = 0; // bytes from one sample to the next
    size_t image_row = 0;     // bytes per row of an image of samples, 0 for a signal (one dimension)
    size_t image_phase = 0;   // position of the first byte of a sample, modulo the stride
    size_t swap = 0;        // byte pairs swapped before coding: 0 none, 1 from byte 0, 2 from byte 1
    bool no_text = false;   // without the text contexts (words, lines), for numeric records
    int level = 0;          // 0: the smallest messages; 1 and 2: faster, a little larger
    int table_bits = 0;     // context tables of 2^table_bits slots
    bool binary = false;    // sparse contexts, for data with control bytes
    int x86 = 0;            // x86 code, branch transform and instruction contexts: 0 none, 1 32-bit, 2 64-bit
};

// Context tables for the data and the level: one bit smaller for the faster levels
int table_size(const size_t data_bytes, const config& cfg, const int forced_table_bits) {
    return forced_table_bits > 0 ? forced_table_bits : table_bits(data_bytes, cfg.level > 0 ? 21 : 22);
}

class model {
public:
    explicit model(const config& cfg): _cfg(cfg), _final(1, 1, 131), _map1(256 * 256, 1311), _map2(size_t{1} << 16, 1311), _map3(size_t{1} << 16, 1311), _match(6) {
        // Byte contexts, left out by the bitmap variant: orders with their counter limits (stable statistics for
        // low orders, adaptive ones above), words, line structure, last occurrences. For text, the faster levels
        // keep fewer orders, in nibble buckets, and leave out the line contexts and the last 3 bytes, which gain
        // least. Folded data spends its time on the fold templates and needs all its orders: the faster levels
        // only leave out its text contexts
        const bool bytes = !cfg.bitmap;
        const bool fast = cfg.level > 0;
        const bool fast_text = fast && !cfg.fold;
        if (bytes) {
            _orders = fast_text ? std::vector<size_t>{0, 1, 2, 3, 4, 6} : std::vector<size_t>{0, 1, 2, 3, 4, 5, 6, 8};
            _limits = fast_text ? std::vector<int>{255, 255, 127, 60, 30, 30} : std::vector<int>{255, 255, 127, 60, 30, 30, 30, 30};
        }
        _text = bytes && !(fast && cfg.fold);
        _full_text = bytes && !fast;
        _buckets = bytes && fast_text;
        _maps_on = cfg.level < 2;
        _record_contexts = bytes && cfg.record > 0 ? 8 : 0;
        if (cfg.sample_bytes > 0) {
            // Pixels: low orders only, no text
            _orders = {0, 1, 2, 3};
            _limits = {255, 255, 127, 60};
            _text = _full_text = false;
            _image_byte_contexts = 6;
            _image_bit_contexts = 11;
        }
        if (cfg.no_text) {
            // Numeric records: no words nor lines
            _text = _full_text = false;
        }
        // Binary data at level 0: sparse contexts. x86 code: instruction contexts too, in place of the text
        // contexts, which gain nothing there once these are in
        _sparse = bytes && !fast && cfg.binary ? (cfg.x86 > 0 ? 5 : 8) : 0;
        _code = bytes && !fast && cfg.x86 > 0 ? 7 : 0;
        if (_code > 0) {
            _text = _full_text = false;
            _parser.emplace(cfg.x86 == 2);
            _last_instruction = std::make_unique<checked_table<last_byte_slot>>(cfg.table_bits - 2);
        }
        _byte_contexts = _orders.size() + _sparse + _code + _record_contexts + (_text ? 2 + 1 + 1 : 0) + (_full_text ? 2 + 1 + 1 : 0) + _image_byte_contexts;

        // Fold contexts: bit templates of three sizes and pixel runs (keyed alone), byte templates with the
        // partial byte, one row above, two rows above, and one layer behind
        if (cfg.fold) {
            const size_t l = cfg.width;
            _bit_templates.push_back({1, 2, l + 1, l, l - 1, 2 * l});
            _bit_templates.push_back({1, 2, 3, 4, l + 2, l + 1, l, l - 1, l - 2, 2 * l + 1, 2 * l, 2 * l - 1, l + 3, l - 3});
            std::vector<size_t> large;
            for (size_t d = 1; d <= 6; ++d) large.push_back(d);
            for (size_t k = 0; k <= 8; ++k) large.push_back(l + 4 - k);
            for (size_t k = 0; k <= 6; ++k) large.push_back(2 * l + 3 - k);
            for (size_t k = 0; k <= 2; ++k) large.push_back(3 * l + 1 - k);
            _bit_templates.push_back(large);
            _byte_templates = bytes ? 2 + (cfg.layer > 0 ? 1 : 0) : 0;

            // Templates read at every bit, compiled once
            for (const auto& distances : _bit_templates) {
                _compiled.emplace_back(distances);
            }
            _runs_template.emplace(std::vector<size_t>{1, l, l - 1});
            _set_template.emplace(std::vector<size_t>{1, 2, l - 1, l, l + 1, 2 * l, l - 2, l + 2});
            // The byte above, two above and one layer behind, aligned on the current byte's first bit: one
            // template per bit position in the byte
            const size_t rows[] = {l, 2 * l, l * cfg.layer};
            for (size_t r = 0; r < _byte_templates; ++r) {
                for (size_t b = 0; b < 8; ++b) {
                    std::vector<size_t> distances;
                    for (size_t k = 0; k < 8; ++k) {
                        distances.push_back(rows[r] + b - k);
                    }
                    _aligned[r].emplace_back(distances);
                }
            }
        }
        _contexts = _byte_contexts + _bit_templates.size() + (cfg.fold ? 1 : 0) + _byte_templates + _image_bit_contexts;

        for (size_t t = 0; t < _contexts; ++t) {
            _tables.emplace_back(cfg.table_bits);
            _maps.emplace_back(1023);
            _limits.resize(std::max(_limits.size(), t + 1), 30);
        }
        _slots.resize(_contexts);
        _hashes.resize(_contexts);
        _nibble_hashes.resize(_contexts);
        if (_text) {
            _last_word = std::make_unique<checked_table<last_byte_slot>>(cfg.table_bits - 2);
        }
        if (_full_text) {
            _last_three = std::make_unique<checked_table<last_byte_slot>>(cfg.table_bits - 2);
        }

        // Mixers with weights selected by: the partial byte, the highest order seen with the match state, the
        // previous byte and the bit position, the previous byte's class (not for text at the faster levels), the
        // fold neighbourhood
        const size_t inputs = _contexts * 3 + 2 + 1;
        _mixers.emplace_back(inputs, 256, 197);
        _mixers.emplace_back(inputs, (_orders.size() + 1) * 4, 197);
        _mixers.emplace_back(inputs, 256 * 8, 197);
        _set_kinds = {0, 1, 2};
        if (!fast_text) {
            _mixers.emplace_back(inputs, 8 * 256, 197);
            _set_kinds.push_back(3);
        }
        if (cfg.fold) {
            _mixers.emplace_back(inputs, 256, 197);
            _set_kinds.push_back(4);
        }
        if (_record_contexts > 0) {
            _mixers.emplace_back(inputs, 256 * 8, 197);
            _set_kinds.push_back(7);
        }
        if (_image_bit_contexts > 0) {
            _mixers.emplace_back(inputs, 33 * 16, 197);
            _set_kinds.push_back(8);
            _mixers.emplace_back(inputs, 33 * 32, 197);
            _set_kinds.push_back(9);
        }
        // Code: weights selected by the instruction field with the partial byte, and by the opcode with the bit
        // position
        if (_code > 0) {
            _mixers.emplace_back(inputs, 64 * 256, 197);
            _set_kinds.push_back(5);
            _mixers.emplace_back(inputs, 1024 * 8, 197);
            _set_kinds.push_back(6);
        }
        _final = logistic_mixer(_mixers.size() + 1, 1, 131, 65536 / static_cast<int>(_mixers.size()));
        _x.reserve(inputs);
        _stream.push_back(0);
        new_byte();
    }

    // Probability that the next bit is 1
    int p() {
        // Mixers learn faster at the start
        const auto rate = 197 + static_cast<int>(uint64_t{197} * 6 * 65536 / (65536 + _bits));
        for (auto& mixer : _mixers) {
            mixer.set_rate(rate);
        }
        _x.clear();

        if (_cfg.fold) {
            fold_contexts();
        }
        if (_image_bit_contexts > 0) {
            sample_contexts();
        }

        // Each context: its probability, its history's probability across contexts, and a deterministic input
        // when only one bit value was ever seen
        const auto input = [&](const int x) {
            _x.push_back(x);
        };
        size_t highest_seen = 0;
        // Nibble buckets: a byte context's slot in a bucket of 16 for each half byte, at the node of the bits
        // coded in that half, instead of a slot anywhere for each bit
        if (_buckets && (_bit_count == 0 || _bit_count == 4)) {
            const uint64_t first = _bit_count == 4 ? _c0 : 0;
            for (size_t t = 0; t < _byte_contexts; ++t) {
                _nibble_hashes[t] = _hashes[t] + first * 0x9E3779B97F4A7C15ull;
            }
        }
        const int half = _bit_count & 3;
        const size_t node = (static_cast<size_t>(_c0) & ((size_t{1} << half) - 1)) | (size_t{1} << half);
        for (size_t t = 0; t < _contexts; ++t) {
            const uint64_t key = t < _byte_contexts ? _hashes[t] ^ (static_cast<uint64_t>(_c0) * 0x9E3779B97F4A7C15ull) : _hashes[t];
            auto& slot = _buckets && t < _byte_contexts ? _tables[t].at(_nibble_hashes[t], node) : _tables[t].at(key);
            _slots[t] = &slot;
            if (!_buckets && t < _byte_contexts && _bit_count < 7) {
                // The slots of the next bit, whichever it is, fetched while this one is computed
                const uint64_t next = static_cast<uint64_t>(_c0) << 1;
                _tables[t].prefetch(_hashes[t] ^ (next * 0x9E3779B97F4A7C15ull));
                _tables[t].prefetch(_hashes[t] ^ ((next | 1) * 0x9E3779B97F4A7C15ull));
            }
            input(slot.seen() ? stretch(clamp_probability(slot.p16())) : 0);
            input(stretch(clamp_probability(_maps[t].p(slot.history))));
            const int zeros = slot.history.zeros(), ones = slot.history.ones();
            input(zeros == 0 && ones > 0 ? ones * 64 : ones == 0 && zeros > 0 ? -zeros * 64 : 0);
            if (slot.seen() && t < _orders.size()) {
                highest_seen = t + 1;
            }
        }

        _expected = _match.expected_bit(_history, _c0, _bit_count);
        const uint32_t length = _match.length();
        _bucket = _expected < 0 ? 0 : std::min(length, 31u) + 1;
        const int sign = _expected == 1 ? 1 : -1;
        input(_expected < 0 ? 0 : sign * stretch(clamp_probability(_match_counters[_bucket].p16())));
        input(_expected < 0 ? 0 : sign * static_cast<int>(std::min(length, 32u)) * 16);
        const size_t match_state = _expected < 0 ? 0 : length < 16 ? 1 : length < 32 ? 2 : 3;
        input(77);

        const size_t sets[] = {
            _c0,
            highest_seen * 4 + match_state,
            static_cast<size_t>(_last) * 8 + static_cast<size_t>(_bit_count),
            static_cast<size_t>(classify(_last)) * 256 + _c0,
            _fold_set,
            _code_state * 256 + _c0,
            _code_opcode * 8 + static_cast<size_t>(_bit_count),
            _column_set * 8 + static_cast<size_t>(_bit_count),
            _image_j * 16 + _activity,
            _image_j * 32 + _image_error,
        };
        // Every mixer reads the same inputs, without copying them
        for (size_t i = 0; i < _mixers.size(); ++i) {
            _final.add(stretch(_mixers[i].mix(_x, sets[_set_kinds[i]])));
        }
        _final.add(77);
        _mixed = _final.mix(0);
        if (!_maps_on) {
            return _mixed;
        }

        // Chained maps: order 1, then order 2 refining the first
        const int p1 = _map1.refine(_mixed, static_cast<size_t>(_last) * 256 + _c0);
        const int p2 = _map2.refine(p1, (hash_mix(context_hash(_history, _history.size(), 2)) ^ _c0) & 0xFFFF);
        if (_image_bit_contexts > 0) {
            // Samples: a third map, on the bit of the sample and where the adaptive prediction stands
            const int p3 = _image_apm.refine(_mixed, _apm_context);
            return clamp_probability((_mixed + p1 + p2 + p3) / 4);
        }
        if (_code == 0) {
            return clamp_probability((_mixed + p1 + 2 * p2) / 4);
        }
        // Code: a third map, on the instruction field and opcode
        const int p3 = _map3.refine(_mixed, (_code_key * 256 + _c0) & 0xFFFF);
        return clamp_probability((_mixed + p1 + p2 + p3) / 4);
    }

    // Learns the bit, after it was coded or decoded
    void update(const int bit) {
        for (auto& mixer : _mixers) {
            mixer.update(bit);
        }
        _final.update(bit);
        for (size_t t = 0; t < _contexts; ++t) {
            _maps[t].update(_slots[t]->history, bit);
            _slots[t]->update(bit, _limits[t]);
        }
        if (_expected >= 0) {
            _match_counters[_bucket].update(bit == _expected, 255);
            if (bit != _expected) {
                _match.miss();
            }
        }
        if (_maps_on) {
            _map1.update(bit);
            _map2.update(bit);
            if (_image_bit_contexts > 0) {
                _image_apm.update(bit);
            }
            if (_code > 0) {
                _map3.update(bit);
            }
        }

        ++_bits;
        _stream.back() |= static_cast<uint8_t>(bit << (7 - _bit_count));
        _c0 = (_c0 << 1) | static_cast<uint32_t>(bit);
        if (++_bit_count == 8) {
            byte_done(static_cast<uint8_t>(_c0 & 0xFF));
        }
    }

private:
    void fold_contexts() {
        const size_t l = _cfg.width;
        const size_t pos = _bits;
        size_t t = _byte_contexts;

        for (size_t k = 0; k < _bit_templates.size(); ++k) {
            _hashes[t++] = _compiled[k].at(_stream, pos) * 0x9E3779B97F4A7C15ull + 977 * (k + 1);
        }
        const auto runs = runs_at(_stream, pos, l);
        const uint64_t runs_key = (std::min<size_t>(runs.left, 63) << 6 | std::min<size_t>(runs.above, 63)) << 3 | _runs_template->at(_stream, pos);
        _hashes[t++] = runs_key * 0x9E3779B97F4A7C15ull + 4444;

        if (_byte_templates > 0) {
            // The byte above (and two above, and one layer behind), aligned on the current byte's first bit
            const auto aligned = [&](const size_t row) {
                return _aligned[row][static_cast<size_t>(_bit_count)].at(_stream, pos);
            };
            const uint64_t above = aligned(0);
            const uint64_t above2 = aligned(1);
            _hashes[t++] = (above << 9 | _c0) * 0x9E3779B97F4A7C15ull + 1111;
            _hashes[t++] = ((above2 << 8 | above) << 9 | _c0) * 0x9E3779B97F4A7C15ull + 2222;
            if (_byte_templates > 2) {
                const uint64_t behind = aligned(2);
                _hashes[t++] = ((behind << 8 | above) << 9 | _c0) * 0x9E3779B97F4A7C15ull + 3333;
            }
        }
        _fold_set = _set_template->at(_stream, pos);
    }

    // Numeric contexts, at every bit of a sample: where each prediction of the sample stands from the middle of
    // the values left possible by the bits known so far, the next bit being 1 when the sample is at least that
    // middle. The statistics are shared by every value: what is learned is the shape of the prediction errors
    void sample_contexts() {
        size_t t = _contexts - _image_bit_contexts;
        const int bits = 8 * static_cast<int>(_cfg.sample_bytes);
        if (!_in_sample) {
            // A byte between samples: one context of its own, the same sets
            _image_j = static_cast<size_t>(bits);
            for (size_t k = 0; k < _image_bit_contexts; ++k) {
                _hashes[t++] = hash_mix((k + 1) << 20 | static_cast<uint64_t>(_bit_count) << 14) + 4444;
            }
            _image_error = 0;
            _apm_context = static_cast<size_t>(bits) * 72;
            return;
        }
        const int j = 8 * static_cast<int>(_sample_byte) + _bit_count;
        const int64_t prefix = _sample_prefix << _bit_count | static_cast<int64_t>(_c0 & ((1u << _bit_count) - 1));
        const int64_t middle = bit_middle(prefix, j, bits);
        _image_j = static_cast<size_t>(j);
        int errors[10];
        for (size_t k = 0; k < 10; ++k) {
            errors[k] = quantize_difference(_predictions[k] - middle);
            _hashes[t++] = hash_mix((k + 1) << 20 | static_cast<uint64_t>(j) << 14 | static_cast<uint64_t>(errors[k] + 128) << 5 | _activity) + 5555;
        }
        _hashes[t++] = hash_mix(uint64_t{9} << 40 | static_cast<uint64_t>(j) << 24 | static_cast<uint64_t>(errors[2] + 128) << 16
                                | static_cast<uint64_t>(errors[0] + 128) << 8 | static_cast<uint64_t>(errors[1] + 128)) + 6666;
        _image_error = static_cast<size_t>(std::max(-15, std::min(16, errors[8] / 2)) + 15);
        _apm_context = static_cast<size_t>(j) * 72 + static_cast<size_t>(std::max(-35, std::min(35, errors[8])) + 36);
    }

    void sample_new_byte(const size_t pos, size_t& t) {
        const size_t stride = _cfg.sample_stride, size = _cfg.sample_bytes;
        const size_t b = (pos + stride - _cfg.image_phase % stride) % stride;
        _in_sample = b < size;
        _sample_byte = b;
        const int bits = 8 * static_cast<int>(size);
        const int64_t top = (int64_t{1} << bits) - 1;
        if (!_in_sample) {
            for (size_t k = 0; k < _image_byte_contexts; ++k) {
                _hashes[t++] = hash_mix(uint64_t{77} << 32 | b << 8 | k) + 100;
            }
            return;
        }
        const size_t start = pos - b;
        _sample_prefix = 0;
        uint64_t known = 0;
        for (size_t i = 0; i < b; ++i) {
            _sample_prefix = _sample_prefix << 8 | _history[start + i];
            known = (known << 8 | _history[start + i]);
        }
        const size_t row = _cfg.image_row;
        const auto sample = [&](const size_t back) -> int64_t {
            if (back == 0 || back > start) {
                return 0;
            }
            int64_t v = 0;
            for (size_t i = 0; i < size; ++i) {
                v = v << 8 | _history[start - back + i];
            }
            return v;
        };
        const auto clip = [&](const int64_t v) { return std::max<int64_t>(0, std::min<int64_t>(top, v)); };
        const auto limit = [](const int64_t v) { return std::max<int64_t>(-(int64_t{1} << 24), std::min<int64_t>(int64_t{1} << 24, v)); };
        const int64_t w = sample(stride), ww = sample(2 * stride), www = sample(3 * stride);
        int64_t activity;
        if (b == 0 && _lms_pending) {
            // The sample just coded, learned before predicting the next one
            const int64_t actual = sample(stride);
            auto* lms = row > 0 ? _lms : _signal_lms;
            lms[0].update(actual - _lms_base[0]);
            lms[1].update(actual - _lms_base[1]);
            _lms_pending = false;
        }
        if (row > 0) {
            const int64_t n = sample(row), nw = sample(row + stride), ne = sample(row - stride);
            const int64_t nn = sample(2 * row), nne = sample(2 * row - stride);
            const int64_t med = std::min(std::max(w, n), std::max(std::min(w, n), w + n - nw));
            _predictions[0] = w;
            _predictions[1] = n;
            _predictions[2] = med;
            _predictions[3] = clip(w + n - nw);
            _predictions[4] = clip(n + ne - nne);
            _predictions[5] = clip(2 * n - nn);
            _predictions[6] = clip(2 * w - ww);
            _predictions[7] = (w + ne + 1) / 2;
            if (b == 0) {
                const int64_t nww = sample(row + 2 * stride), nnw = sample(2 * row + stride), nnn = sample(3 * row), nee = sample(row - 2 * stride);
                const int64_t f0[] = {n - w, nw - w, ne - w, ww - w, nn - w, nne - w, nww - w, nnw - w, nee - w, nnn - w};
                const int64_t f1[] = {w - n, nw - n, ne - n, ww - n, nn - n, nne - n, nww - n, nnw - n, nee - n, nnn - n};
                int64_t g0[10], g1[10];
                for (size_t i = 0; i < 10; ++i) {
                    g0[i] = limit(f0[i]);
                    g1[i] = limit(f1[i]);
                }
                _lms_base[0] = w;
                _lms_base[1] = n;
                _lms_value[0] = clip(w + _lms[0].predict(g0));
                _lms_value[1] = clip(n + _lms[1].predict(g1));
                _lms_pending = start >= 3 * row + 2 * stride;
            }
            activity = std::abs(w - nw) + std::abs(n - nw) + std::abs(n - ne) + std::abs(w - ww) / 2 + std::abs(n - nn) / 2;
        } else {
            // A signal: extrapolations of the last samples
            const int64_t w4 = sample(4 * stride);
            _predictions[0] = w;
            _predictions[1] = ww;
            _predictions[2] = clip(2 * w - ww);
            _predictions[3] = clip(3 * w - 3 * ww + www);
            _predictions[4] = (w + ww + 1) / 2;
            _predictions[5] = clip(w + (w - ww) / 2);
            _predictions[6] = www;
            _predictions[7] = (w + ww + www + w4 + 2) / 4;
            if (b == 0) {
                const size_t taps = _signal_lms[0].size();
                std::vector<int64_t>& f = _signal_features;
                f.resize(taps);
                for (size_t i = 0; i < taps; ++i) {
                    f[i] = limit(sample((i + 2) * stride) - w);
                }
                _lms_base[0] = _lms_base[1] = w;
                _lms_value[0] = clip(w + _signal_lms[0].predict(f.data()));
                _lms_value[1] = clip(w + _signal_lms[1].predict(f.data()));
                _lms_pending = start >= (taps + 2) * stride;
            }
            activity = std::abs(w - ww) + std::abs(ww - www);
        }
        _predictions[8] = _lms_value[0];
        _predictions[9] = _lms_value[1];
        const int scale = std::max(0, bits - 16);   // activity on the scale of 16-bit samples
        _activity = static_cast<size_t>(std::min(15, 64 - (activity >> scale == 0 ? 64 : leading_zeros(static_cast<uint64_t>(activity >> scale)))));
        const uint64_t h = hash_mix((b + 1) << 8 ^ known * 0x9E3779B97F4A7C15ull);
        const int shift = std::max(0, bits - 16);
        _hashes[t++] = hash_mix(h + static_cast<uint64_t>(w >> shift)) + 101;
        _hashes[t++] = hash_mix(h + static_cast<uint64_t>(_predictions[1] >> shift)) + 102;
        _hashes[t++] = hash_mix(h + static_cast<uint64_t>(_predictions[2] >> shift)) + 103;
        _hashes[t++] = hash_mix(h + static_cast<uint64_t>(_predictions[3] >> shift)) + 104;
        _hashes[t++] = hash_mix(h + (static_cast<uint64_t>(_predictions[7] >> (shift + 2)) << 4 | _activity)) + 105;
        _hashes[t++] = hash_mix(h + static_cast<uint64_t>(_predictions[8] >> shift)) + 106;
    }

    void byte_done(const uint8_t byte) {
        _stream.push_back(0);
        _history.push_back(byte);
        _last = byte;
        _lines.byte_added(_history);
        // Capitals escapes neither extend nor end a word
        if (byte != 0x01 && byte != 0x02) {
            if (std::isalpha(byte)) {
                _word = hash_mix(_word ^ static_cast<uint64_t>(std::tolower(byte)) * 0x2545F4914F6CDD1Dull);
            } else if (_word != 0) {
                _previous_word = _word;
                _word = 0;
            }
        }
        _match.byte_added(_history);
        if (_code > 0) {
            _parser->byte_added(byte);
            _instruction_slot->update(byte);
            _instruction_slot2->update(byte);
        }
        if (_text) {
            _word_slot->update(byte);
        }
        if (_full_text) {
            _three_slot->update(byte);
        }
        new_byte();
    }

    // Binary data: bytes 2 and 3 back alone, the high nibbles of the last 4, bytes 4 and 8 back, the 2 last bytes
    // with the position in a 32-bit word; and for data that is not x86 code, contexts of 32-bit words: the
    // bytes of the current word so far with the byte 4 back, or with the whole previous word, and the bytes 4
    // and 8 back with the position in the word
    void sparse_contexts(size_t& t, const size_t pos) {
        const auto& h = _history;
        const uint64_t column = pos & 3;
        _hashes[t++] = sparse_hash(h, pos, {2});
        _hashes[t++] = sparse_hash(h, pos, {3});
        _hashes[t++] = sparse_hash(h, pos, {1, 2, 3, 4}, 0xF0);
        _hashes[t++] = sparse_hash(h, pos, {4, 8});
        _hashes[t++] = sparse_hash(h, pos, {1, 2}) + column;
        if (_sparse > 5) {
            static const std::initializer_list<size_t> word[4] = {{4}, {1, 4}, {1, 2, 4}, {1, 2, 3, 4}};
            static const std::initializer_list<size_t> words[4] = {{1, 2, 3, 4}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6, 7}};
            _hashes[t++] = sparse_hash(h, pos, word[column]) + column * 11;
            _hashes[t++] = sparse_hash(h, pos, words[column]) + column * 13;
            _hashes[t++] = sparse_hash(h, pos, {4, 8}) + column * 17 + 1;
        }
    }

    // x86 code: the field of the next byte in its instruction (opcode, ModRM, SIB, displacement, immediate) and
    // its index there, with the opcode (and the prefixes, or the last ones while at the opcode), the ModRM once
    // past it, the previous instructions' opcodes and ModRM, and what followed the last time in the same place
    void code_contexts(size_t& t, const size_t pos) {
        const auto& x = *_parser;
        const uint64_t field = x.next_field(), index = static_cast<uint64_t>(std::min(x.index(), 7));
        const uint64_t state = field * 8 + index, op = x.current();
        const bool at_opcode = field == x86_parser::opcode;
        const uint64_t at = at_opcode ? x.prefixes() : op | uint64_t{x.prefixes() & 0xFF} << 12;
        const uint64_t modrm = field >= x86_parser::sib ? x.modrm_byte() : 0;
        const uint64_t p1 = x.previous(1), p2 = x.previous(2), p3 = x.previous(3), pm = x.previous_modrm();
        const uint64_t b1 = pos >= 1 ? _history[pos - 1] : 0, b2 = pos >= 2 ? _history[pos - 2] : 0;
        const uint64_t place = state << 24 | at << 8 | modrm;

        const uint64_t key = place * 0x9E3779B97F4A7C15ull;
        _instruction_slot = &_last_instruction->at(key + 0x31);
        _instruction_slot2 = &_last_instruction->at(key + b1 * 0x2545F4914F6CDD1Dull + 0x77);
        const uint64_t last1 = _instruction_slot->count > 0 ? _instruction_slot->byte : 256u;
        const uint64_t last2 = _instruction_slot2->count > 0 ? _instruction_slot2->byte : 256u;

        const uint64_t hashes[] = {
            place,
            at_opcode ? (state << 40 | p1 << 20 | p2 | uint64_t{1} << 60) : place << 8 | b1,
            state << 40 | at << 20 | p1,
            at_opcode ? (p1 << 40 | p2 << 20 | p3 | uint64_t{1} << 60) : place << 16 | b1 << 8 | b2,
            state << 40 | at << 20 | p1 << 8 | pm,
            place << 9 | last1,
            (state << 9 | last1) << 9 | last2,
        };
        for (size_t k = 0; k < 7; ++k) {
            _hashes[t++] = hash_mix(hashes[k] + (k + 1) * 0x51ED27) + 4242;
        }
        _code_state = static_cast<size_t>(std::min<uint64_t>(state, 63));
        _code_opcode = static_cast<size_t>(op & 0x3FF);
        _code_key = hash_mix(state << 16 | (at_opcode ? p1 : op));
    }

    void new_byte() {
        _c0 = 1;
        _bit_count = 0;
        const size_t pos = _history.size();
        size_t t = 0;
        for (const auto order : _orders) {
            _hashes[t++] = context_hash(_history, pos, order);
        }
        if (_image_byte_contexts > 0) {
            sample_new_byte(pos, t);
        }
        if (_record_contexts > 0) {
            const size_t r = _cfg.record;
            const uint64_t c = pos % r;
            const auto at = [&](const size_t back) -> uint64_t {
                return pos >= back ? _history[pos - back] : 256;
            };
            const uint64_t above = at(r), above2 = at(2 * r), prev = at(1), prev2 = at(2);
            const uint64_t ne = at(r - 1), nw = at(r + 1);
            _column_set = std::min<size_t>(c, 255);
            const uint64_t base = (c + 1) * 0x2545F4914F6CDD1Dull;
            _hashes[t++] = hash_mix(base + above) + 1;
            _hashes[t++] = hash_mix(base + (above << 9 | above2)) + 2;
            _hashes[t++] = hash_mix(base + prev) + 3;
            _hashes[t++] = hash_mix(base + (prev << 9 | above)) + 4;
            _hashes[t++] = hash_mix(base + (above << 18 | ne << 9 | nw)) + 5;
            _hashes[t++] = hash_mix(base + ((2 * above - above2) & 0x3FF)) + 6;
            _hashes[t++] = hash_mix(base) + 7;
            _hashes[t++] = hash_mix(base + (prev << 9 | prev2)) + 8;
        }
        if (_sparse > 0) {
            sparse_contexts(t, pos);
        }
        if (_code > 0) {
            code_contexts(t, pos);
        }
        if (!_text) {
            return;
        }

        // Words: the current one, the previous one with it
        _hashes[t++] = _word * 31 + 101;
        _hashes[t++] = hash_mix(_previous_word * 17 + _word) + 202;

        // Line structure: the byte above, the column, the previous line's length
        const uint64_t last = pos > 0 ? _history[pos - 1] : 256;
        const size_t column = _lines.column(_history);
        const size_t previous = _lines.previous_length();
        if (_full_text) {
            const uint64_t above = _lines.above(_history);
            _hashes[t++] = hash_mix(above * 257 + last) + 707;
            _hashes[t++] = hash_mix(above * 64 + std::min<size_t>(column, 63)) + 808;
        }
        _hashes[t++] = hash_mix(std::min<size_t>(column, 80) * 1031 + std::min<size_t>(previous, 120) * 3 + (last == ' ')) + 909;
        if (_full_text) {
            const long to_end = static_cast<long>(column) - static_cast<long>(previous);
            _hashes[t++] = hash_mix(static_cast<uint64_t>(std::max(-40L, std::min(40L, to_end)) + 40) * 8
                                    + static_cast<uint64_t>(classify(static_cast<uint8_t>(last)))) + 1010;
        }

        // What followed the last occurrence of the current word, and of the last 3 bytes
        _word_slot = &_last_word->at(_word + 0x77);
        _hashes[t++] = hash_mix((_word_slot->count > 0 ? _word_slot->byte : 256u) * 257 + last) + 1212;
        if (_full_text) {
            _three_slot = &_last_three->at(context_hash(_history, pos, 3) + 0x99);
            _hashes[t++] = hash_mix((_three_slot->count > 0 ? _three_slot->byte : 256u) * 257 + last) + 1313;
        }
    }

    const config& _cfg;
    std::vector<size_t> _orders;
    std::vector<int> _limits;
    size_t _byte_contexts = 0, _contexts = 0, _byte_templates = 0, _record_contexts = 0, _column_set = 0;
    size_t _image_byte_contexts = 0, _image_bit_contexts = 0, _image_j = 0, _activity = 0, _image_error = 0;
    bool _in_sample = false;
    size_t _sample_byte = 0;
    int64_t _sample_prefix = 0;
    int64_t _predictions[10] = {};
    lms_predictor _lms[2] = {lms_predictor(10, 3), lms_predictor(10, 6)};
    int64_t _lms_base[2] = {}, _lms_value[2] = {};
    lms_predictor _signal_lms[2] = {lms_predictor(16, 1), lms_predictor(16, 4)};
    std::vector<int64_t> _signal_features;
    bool _lms_pending = false;
    size_t _apm_context = 0;
    probability_map _image_apm{33 * 72, 1311};
    bool _text = false, _full_text = false, _buckets = false, _maps_on = true;
    std::vector<size_t> _set_kinds;
    std::vector<std::vector<size_t>> _bit_templates;
    std::vector<compiled_template> _compiled, _aligned[3];
    std::optional<compiled_template> _runs_template, _set_template;
    std::vector<checked_table<context_slot>> _tables;
    std::vector<history_map> _maps;
    std::vector<context_slot*> _slots;
    std::vector<uint64_t> _hashes, _nibble_hashes;
    std::unique_ptr<checked_table<last_byte_slot>> _last_word, _last_three;
    last_byte_slot *_word_slot = nullptr, *_three_slot = nullptr;
    std::vector<logistic_mixer> _mixers;
    logistic_mixer _final;
    std::vector<int> _x;
    probability_map _map1, _map2, _map3;
    match_model _match;
    bit_counter _match_counters[33];
    line_tracker _lines;
    std::vector<uint8_t> _history, _stream;
    size_t _bits = 0, _bucket = 0, _fold_set = 0;
    uint32_t _c0 = 1;
    int _bit_count = 0, _expected = -1;
    uint8_t _last = 0;
    uint64_t _word = 0, _previous_word = 0;
    size_t _sparse = 0, _code = 0, _code_state = 0, _code_opcode = 0;
    uint64_t _code_key = 0;
    std::optional<x86_parser> _parser;
    std::unique_ptr<checked_table<last_byte_slot>> _last_instruction;
    last_byte_slot *_instruction_slot = nullptr, *_instruction_slot2 = nullptr;
    int _mixed = 32768;
};

void encode_number(arithmetic_encoder& encoder, const size_t value, const int bits) {
    for (int i = bits - 1; i >= 0; --i) {
        encoder.encode(static_cast<int>((value >> i) & 1), 32768);
    }
}

size_t decode_number(arithmetic_decoder& decoder, const int bits) {
    size_t value = 0;
    for (int i = 0; i < bits; ++i) {
        value = (value << 1) | static_cast<size_t>(decoder.decode(32768));
    }
    return value;
}

// Header, at p = 1/2: payload size (32 bits), fold (1), bitmap (1), capitals (1), fold width (16), layer (8),
// level (2), binary (1), x86 (2). Then
// the payload, and the CRC-32 of the original data (32 bits), so that decoding to anything else is detected
gmprime compress(const std::vector<uint8_t>& payload, config cfg, const bool capitals, const int forced_table_bits,
                 const uint32_t check) {
    arithmetic_encoder encoder;
    encode_number(encoder, payload.size(), 32);
    encode_number(encoder, cfg.fold ? 1 : 0, 1);
    encode_number(encoder, cfg.bitmap ? 1 : 0, 1);
    encode_number(encoder, capitals ? 1 : 0, 1);
    encode_number(encoder, cfg.width, 16);
    encode_number(encoder, cfg.layer, 8);
    encode_number(encoder, static_cast<size_t>(cfg.level), 2);
    encode_number(encoder, cfg.binary ? 1 : 0, 1);
    encode_number(encoder, static_cast<size_t>(cfg.x86), 2);
    // The structure found beyond the fold, only when there is one
    const bool structured = cfg.width >= 65536 || cfg.record > 0 || cfg.lattice_width > 0 || cfg.sample_bytes > 0 || cfg.swap > 0 || cfg.no_text;
    encode_number(encoder, structured ? 1 : 0, 1);
    if (structured) {
        encode_number(encoder, cfg.width >> 16, 8);
        encode_number(encoder, cfg.record, 16);
        encode_number(encoder, cfg.lattice_width, 16);
        encode_number(encoder, cfg.sample_bytes, 4);
        encode_number(encoder, cfg.sample_stride, 16);
        encode_number(encoder, cfg.image_row, 24);
        encode_number(encoder, cfg.image_phase, 16);
        encode_number(encoder, cfg.swap, 2);
        encode_number(encoder, cfg.no_text ? 1 : 0, 1);
        encode_number(encoder, cfg.lattice.size(), 4);
        for (const auto& f : cfg.lattice) {
            encode_number(encoder, f.offset, 16);
            encode_number(encoder, f.format.size == 8 ? 1 : 0, 1);
            encode_number(encoder, f.format.big_endian ? 1 : 0, 1);
            encode_number(encoder, f.format.ibm ? 1 : 0, 1);
            encode_number(encoder, f.offset_binary ? 1 : 0, 1);
            encode_number(encoder, f.step_bits >> 32, 32);
            encode_number(encoder, f.step_bits & 0xFFFFFFFFu, 32);
            encode_number(encoder, static_cast<size_t>(f.residual_bits), 3);
            encode_number(encoder, f.delta ? 1 : 0, 1);
        }
    }

    cfg.table_bits = table_size(payload.size(), cfg, forced_table_bits);
    model m(cfg);
    for (const auto byte : payload) {
        for (int i = 7; i >= 0; --i) {
            const int bit = (byte >> i) & 1;
            encoder.encode(bit, m.p());
            m.update(bit);
        }
    }
    encode_number(encoder, check, 32);
    return encoder.finish();
}

// The original data, or nothing if it does not match its checksum: a corrupted message, or decompressed with
// another table size than it was compressed with
std::optional<std::vector<uint8_t>> decompress(const gmprime& code, const int forced_table_bits) {
    arithmetic_decoder decoder(code);
    config cfg;
    const size_t size = decode_number(decoder, 32);
    cfg.fold = decode_number(decoder, 1) == 1;
    cfg.bitmap = decode_number(decoder, 1) == 1;
    const bool capitals = decode_number(decoder, 1) == 1;
    cfg.width = decode_number(decoder, 16);
    cfg.layer = decode_number(decoder, 8);
    cfg.level = static_cast<int>(decode_number(decoder, 2));
    cfg.binary = decode_number(decoder, 1) == 1;
    cfg.x86 = static_cast<int>(decode_number(decoder, 2));
    if (decode_number(decoder, 1) == 1) {
        cfg.width |= decode_number(decoder, 8) << 16;
        cfg.record = decode_number(decoder, 16);
        cfg.lattice_width = decode_number(decoder, 16);
        cfg.sample_bytes = decode_number(decoder, 4);
        cfg.sample_stride = decode_number(decoder, 16);
        cfg.image_row = decode_number(decoder, 24);
        cfg.image_phase = decode_number(decoder, 16);
        cfg.swap = decode_number(decoder, 2);
        cfg.no_text = decode_number(decoder, 1) == 1;
        cfg.lattice.resize(decode_number(decoder, 4));
        for (auto& f : cfg.lattice) {
            f.offset = decode_number(decoder, 16);
            f.format.size = decode_number(decoder, 1) ? 8 : 4;
            f.format.big_endian = decode_number(decoder, 1) == 1;
            f.format.ibm = decode_number(decoder, 1) == 1;
            f.offset_binary = decode_number(decoder, 1) == 1;
            f.step_bits = decode_number(decoder, 32) << 32;
            f.step_bits |= decode_number(decoder, 32);
            f.residual_bits = static_cast<int>(decode_number(decoder, 3));
            f.delta = decode_number(decoder, 1) == 1;
        }
    }

    cfg.table_bits = table_size(size, cfg, forced_table_bits);
    model m(cfg);
    std::vector<uint8_t> payload;
    payload.reserve(size);
    for (size_t n = 0; n < size; ++n) {
        int byte = 0;
        for (int i = 0; i < 8; ++i) {
            const int bit = decoder.decode(m.p());
            m.update(bit);
            byte = (byte << 1) | bit;
        }
        payload.push_back(static_cast<uint8_t>(byte));
    }
    const auto check = static_cast<uint32_t>(decode_number(decoder, 32));
    auto data = capitals ? capitals_decode(payload) : std::move(payload);
    if (cfg.x86 > 0) {
        x86_branches(data, false);
    }
    if (cfg.swap > 0) {
        swap_pairs(data, cfg.swap - 1);
    }
    if (cfg.lattice_width > 0) {
        auto original = lattice_decode(data, cfg.lattice_width, cfg.lattice);
        if (!original) {
            return std::nullopt;
        }
        data = std::move(*original);
    }
    if (crc32(data) != check) {
        return std::nullopt;
    }
    return data;
}

struct choice {
    gmprime code;
    config cfg;
    bool capitals = false;
};

// The ways to code the data that the analyses suggest, each with the transforms it needs
struct plan {
    const char* name;
    config cfg;
    bool capitals = false;
};

std::vector<uint8_t> transformed(const std::vector<uint8_t>& data, const plan& p) {
    if (p.capitals) {
        return *capitals_encode(data);
    }
    auto result = p.cfg.lattice_width > 0 ? lattice_encode(data, p.cfg.lattice_width, p.cfg.lattice) : data;
    if (p.cfg.swap > 0) {
        swap_pairs(result, p.cfg.swap - 1);
    }
    return result;
}

// The plans for data that is not x86 code: text (with the capitals transform) unless the fold is strong, the fold
// of analyze_folds (and its bitmap variant when strong), records of the width where bytes repeat with column contexts, their floating point
// fields on a lattice (and each 4-byte one as a signal), 16-bit samples of an image
std::vector<plan> structure_plans(const gmprime& file, const std::vector<uint8_t>& data, const config& base, const std::string& capitals_mode) {
    std::vector<plan> plans;
    // The first 128 KB find the same folds as the whole data, much sooner. A strong fold is not text (the first
    // bytes can say otherwise, on nci, where the fold gains later)
    const auto folds = file.analyze_folds(64, 8192, 131072);
    if (!folds.strong()) {
        plans.push_back({"text", base});
        if (capitals_mode != "off" && base.level == 0 && capitals_encode(data).has_value()) {
            plans.push_back({"capitals", base, true});
        }
    }
    if (folds.width > 0) {
        config c = base;
        c.fold = true;
        c.width = folds.width;
        c.layer = folds.layer < 256 ? folds.layer : 0;
        plans.push_back({"fold", c});
        if (folds.strong()) {
            c.bitmap = true;
            plans.push_back({"bitmap", c});
        }
    }

    const auto records = find_record_width(data, 2048, 131072);
    if (records.width >= 2) {
        config c = base;
        c.fold = true;
        c.width = 8 * records.width;
        c.record = records.width;
        plans.push_back({"records", c});
        const auto fields = find_lattice_fields(data, records.width);
        if (!fields.empty()) {
            c.lattice_width = records.width;
            c.lattice = fields;
            c.record = records.width + 1;
            c.width = 8 * c.record;
            c.no_text = true;
            plans.push_back({"lattice", c});
            // A 4-byte field as a signal: k in offset binary, predicted from the same field in the records before
            for (size_t j = 0; j < fields.size(); ++j) {
                if (fields[j].format.size == 4) {
                    config d = c;
                    d.lattice[j].delta = false;
                    d.lattice[j].offset_binary = true;
                    d.sample_bytes = 4;
                    d.sample_stride = d.record;
                    d.image_phase = 1 + d.lattice[j].offset;
                    plans.push_back({"signal", d});
                }
            }
        }
    }

    // Rows of up to 8192 samples, or an eighth of the data for small images
    const auto image = find_image16(data, 16, std::min<size_t>(8192, data.size() / 8));
    if (image.found && image.row_bytes > 0) {
        config c = base;
        c.fold = true;
        c.width = 8 * image.row_bytes;
        c.image_row = image.row_bytes;
        c.sample_bytes = 2;
        c.sample_stride = 2;
        c.swap = image.little_endian ? image.start + 1 : 0;
        c.image_phase = image.start;
        plans.push_back({"image16", c});
    }
    return plans;
}

// Each plan compresses the first sample bytes, the smallest message wins, then the whole data is compressed that
// way (data that fits in the sample: the smallest message is kept). The folded plans take about twice the time
// of the text ones: they must be 1% smaller to win over them
choice compress_plans(const std::vector<uint8_t>& data, const std::vector<plan>& plans, const int forced_table_bits, const size_t sample,
                      const bool verbose) {
    const uint32_t check = crc32(data);
    const bool whole = sample == 0 || data.size() <= sample;
    const std::vector<uint8_t> first = whole ? data : std::vector<uint8_t>(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(sample));
    choice best;
    size_t best_size = 0;
    for (const auto& p : plans) {
        auto code = compress(transformed(first, p), p.cfg, p.capitals, forced_table_bits, whole ? check : crc32(first));
        const size_t size = code.bytes().size();
        if (verbose) {
            std::fprintf(stderr, "  %-8s %zu bytes\n", p.name, size);
        }
        // Weighed: a folded plan against a text one counts 1% larger
        const bool slower = p.cfg.fold && best_size > 0 && !best.cfg.fold;
        if (best_size == 0 || (slower ? size + size / 100 < best_size : size < best_size)) {
            best_size = size;
            best.cfg = p.cfg;
            best.capitals = p.capitals;
            if (whole) {
                best.code = std::move(code);
            }
        }
    }
    if (!whole) {
        best.code = compress(transformed(data, {"", best.cfg, best.capitals}), best.cfg, best.capitals, forced_table_bits, check);
    }
    return best;
}

// Analyzes the data, then compresses it the way that gives the smaller message: with or without the capitals
// transform for text, with the full or the bitmap model for folded data, unless the options impose one. The
// two ways are compared on the first sample bytes, then the whole data is compressed once the better way;
// with sample 0, or data that fits in the sample, both are compressed whole and the smaller message kept.
// In auto mode, binary data that is not x86 code, and data that folds, is compressed the way of the best structure
// plan instead.
choice compress_best(const gmprime& file, const std::string& mode, const std::string& capitals_mode, const int forced_table_bits,
                     const size_t sample, const int level, const std::string& code_mode, const bool verbose) {
    auto data = file.bytes();
    const uint32_t check = crc32(data);
    choice best;
    best.cfg.level = level;

    // Binary data: 10% of control bytes (below, mostly text, the sparse contexts cost more time than they gain). x86 code: 0.5% of CALL or JMP opcodes with a near target, 64-bit if
    // REX.W prefixes are common; its branches are made absolute
    best.cfg.binary = code_mode == "auto" ? control_density(data) >= 1000 : code_mode != "off";
    best.cfg.x86 = code_mode == "32" ? 1 : code_mode == "64" ? 2 : code_mode != "auto" ? 0
                 : best.cfg.binary && x86_branch_density(data) >= 50 ? (x86_long_mode(data) ? 2 : 1) : 0;
    // Binary data, or data that folds, may hide records or images: each way the analysis suggests is tried. Text
    // keeps the way below, which only compares the capitals transform
    if (mode == "auto" && best.cfg.x86 == 0 && capitals_mode != "on" && (best.cfg.binary || file.analyze_folds(64, 8192, 131072).strong())) {
        return compress_plans(data, structure_plans(file, data, best.cfg, capitals_mode), forced_table_bits, sample, verbose);
    }
    if (best.cfg.x86 > 0) {
        x86_branches(data, true);
    }

    // The first 128 KB find the same folds as the whole data, much sooner
    const auto folds = file.analyze_folds(64, 8192, 131072);
    best.cfg.fold = mode == "fold" || mode == "fold-full" || mode == "bitmap" || (mode == "auto" && folds.strong());
    if (best.cfg.fold) {
        best.cfg.width = folds.width;
        best.cfg.layer = folds.layer < 256 ? folds.layer : 0;
    }
    best.cfg.bitmap = mode == "bitmap";

    const auto lowered = best.cfg.x86 > 0 ? std::nullopt : capitals_encode(data);
    best.capitals = capitals_mode == "on" && !best.cfg.fold && lowered.has_value();
    // The faster levels skip comparing capitals, which compresses the sample twice
    const bool compare_capitals = capitals_mode == "auto" && level == 0 && !best.cfg.fold && lowered.has_value();
    const bool compare_bitmap = best.cfg.fold && (mode == "auto" || mode == "fold");

    if (compare_capitals || compare_bitmap) {
        const bool whole = sample == 0 || data.size() <= sample;
        const auto first = [&](const std::vector<uint8_t>& bytes) {
            return whole ? bytes : std::vector<uint8_t>(bytes.begin(), bytes.begin() + static_cast<std::ptrdiff_t>(sample));
        };
        config variant = best.cfg;
        variant.bitmap = compare_bitmap;
        auto plain = compress(first(data), best.cfg, false, forced_table_bits, check);
        // The capitals transform of the same part of the data: its first sample bytes would hold less text
        auto other = compare_capitals ? compress(*capitals_encode(first(data)), best.cfg, true, forced_table_bits, check)
                                      : compress(first(data), variant, false, forced_table_bits, check);

        if (other.bytes().size() < plain.bytes().size()) {
            best.capitals = compare_capitals;
            best.cfg = variant;
            plain = std::move(other);
        }
        if (whole) {
            best.code = std::move(plain);
            return best;
        }
    }

    best.code = compress(best.capitals ? *lowered : data, best.cfg, best.capitals, forced_table_bits, check);
    return best;
}

int usage() {
    std::fprintf(stderr, "usage: gmz [-1|-2] [-m auto|text|fold|fold-full|bitmap] [-C auto|on|off] [-x auto|off|bin|32|64] [-S bytes]\n"
                         "           [-T bits] c|d INPUT OUTPUT\n"
                         "       gmz [options] t INPUT\n");
    return 2;
}

} // namespace

int main(int argc, char** argv) {
    std::string mode = "auto", capitals_mode = "auto", code_mode = "auto";
    int forced_table_bits = 0;
    size_t sample = 65536;
    int level = 0;
    int arg = 1;
    while (arg + 1 < argc && argv[arg][0] == '-') {
        const std::string flag = argv[arg];
        if (flag == "-1" || flag == "-2") {
            level = flag[1] - '0';
            ++arg;
            continue;
        }
        const char* value = argv[arg + 1];
        if (flag == "-m") mode = value;
        else if (flag == "-C") capitals_mode = value;
        else if (flag == "-x" && (std::string(value) == "auto" || std::string(value) == "off" || std::string(value) == "bin"
                                  || std::string(value) == "32" || std::string(value) == "64")) code_mode = value;
        else if (flag == "-S") sample = std::stoul(value);
        else if (flag == "-T") forced_table_bits = std::stoi(value);
        else return usage();
        arg += 2;
    }
    if (arg >= argc) {
        return usage();
    }
    const std::string command = argv[arg];
    const bool writes = command == "c" || command == "d";
    if ((writes && argc - arg != 3) || (command == "t" && argc - arg != 2) || (!writes && command != "t")) {
        return usage();
    }

    gmprime input;
    input.load_from_file(argv[arg + 1]);

    if (command == "d") {
        const auto data = decompress(input, forced_table_bits);
        if (!data) {
            std::fprintf(stderr, "gmz: the decoded data does not match its checksum\n");
            return 1;
        }
        gmprime::from_bytes(*data).save_to_file(argv[arg + 2]);
        return 0;
    }

    const auto start = std::chrono::steady_clock::now();
    // The test reports the sizes of the plans compared, on standard error
    const auto best = compress_best(input, mode, capitals_mode, forced_table_bits, sample, level, code_mode, command == "t");
    if (command == "c") {
        best.code.save_to_file(argv[arg + 2]);
        return 0;
    }

    // Test: decode what was coded and report
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const auto data = input.bytes();
    const bool same = decompress(best.code, forced_table_bits) == std::optional(data);
    const size_t size = best.code.bytes().size();
    std::printf("%zu -> %zu bytes, %.3f bits/byte, %s, %.1f s, %s%s", data.size(), size, data.empty() ? 0.0 : 8.0 * size / data.size(),
                same ? "round trip ok" : "ROUND TRIP FAILED", seconds, best.cfg.fold ? "fold" : "text", best.capitals ? " + capitals" : "");
    std::printf("%s%s", best.cfg.binary ? " + binary" : "", best.cfg.x86 == 1 ? " + x86-32" : best.cfg.x86 == 2 ? " + x86-64" : "");
    if (best.cfg.fold) {
        std::printf(" (width %zu bits, layer %zu%s)", best.cfg.width, best.cfg.layer, best.cfg.bitmap ? ", bitmap" : "");
    }
    if (best.cfg.record > 0) {
        std::printf(" records %zu", best.cfg.record);
    }
    if (best.cfg.lattice_width > 0) {
        std::printf(" lattice %zu fields", best.cfg.lattice.size());
    }
    if (best.cfg.sample_bytes > 0) {
        std::printf(" samples %zu bytes every %zu, row %zu%s", best.cfg.sample_bytes, best.cfg.sample_stride, best.cfg.image_row, best.cfg.swap ? ", swapped" : "");
    }
    std::printf("\n");
    return same ? 0 : 1;
}
