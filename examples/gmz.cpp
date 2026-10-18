// gmz: a context-mixing compressor built on gmprime and its coding toolkit.
//
//   gmz [options] c INPUT OUTPUT    compress
//   gmz [options] d INPUT OUTPUT    decompress
//   gmz [options] t INPUT           compress, decompress, check and report, without writing anything
//
//   -m MODE   auto (default), text, fold (full and bitmap variants tried), fold-full, bitmap
//   -C MODE   capitals transform for text: auto (tried both ways, default), on, off
//   -S BYTES  ways compared on the first BYTES bytes, then the data compressed once (default 65536); 0 compresses
//             the data both ways and keeps the smaller message, slower and sometimes a little smaller
//   -T BITS   context table size in bits of slot count, instead of one sized for the data (the decompressor
//             needs the same value)
//   -1, -2    faster levels: text compresses 2 to 2.5 times faster and decompresses 1.6 to 1.9 times faster, for
//             messages 1.1% and 1.8% larger over the Canterbury corpus. Fewer contexts and mixers, nibble buckets,
//             smaller tables, no capitals comparison; folded data only loses its text contexts; -2 also drops the
//             chained maps. The level is stored in the message for the decompressor
//
// The file is first analyzed with gmprime::analyze_folds. Data that folds (bitmap rows, fixed-size records)
// is coded with bit and byte templates across the fold, in 2D and 3D; anything else with the text model:
// context orders, words, line structure, last occurrences, indirect modeling, four mixers and chained maps.
// The compressed message is a gmprime: the arithmetic coder's output, saved and loaded like any number.
#include <chrono>
#include <cstdio>
#include <memory>
#include <optional>
#include <string>
#include <gmprime_coding.hpp>

using namespace GMPrime;
using namespace GMPrime::coding;

namespace {

// What the model uses, decided by the analysis and stored in the message header
struct config {
    bool fold = false;      // folded model: templates across the fold
    bool bitmap = false;    // folded model without byte contexts, for pixels
    size_t width = 0;       // fold width in bits
    size_t layer = 0;       // rows per layer of the 3D view, 0 for none
    int level = 0;          // 0: the smallest messages; 1 and 2: faster, a little larger
    int table_bits = 0;     // context tables of 2^table_bits slots
};

// Context tables for the data and the level: one bit smaller for the faster levels
int table_size(const size_t data_bytes, const config& cfg, const int forced_table_bits) {
    return forced_table_bits > 0 ? forced_table_bits : table_bits(data_bytes, cfg.level > 0 ? 21 : 22);
}

class model {
public:
    explicit model(const config& cfg): _cfg(cfg), _final(1, 1, 131), _map1(256 * 256, 1311), _map2(size_t{1} << 16, 1311), _match(6) {
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
        _byte_contexts = _orders.size() + (_text ? 2 + 1 + 1 : 0) + (_full_text ? 2 + 1 + 1 : 0);

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
        _contexts = _byte_contexts + _bit_templates.size() + (cfg.fold ? 1 : 0) + _byte_templates;

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
        return clamp_probability((_mixed + p1 + 2 * p2) / 4);
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
        if (_text) {
            _word_slot->update(byte);
        }
        if (_full_text) {
            _three_slot->update(byte);
        }
        new_byte();
    }

    void new_byte() {
        _c0 = 1;
        _bit_count = 0;
        const size_t pos = _history.size();
        size_t t = 0;
        for (const auto order : _orders) {
            _hashes[t++] = context_hash(_history, pos, order);
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
    size_t _byte_contexts = 0, _contexts = 0, _byte_templates = 0;
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
    probability_map _map1, _map2;
    match_model _match;
    bit_counter _match_counters[33];
    line_tracker _lines;
    std::vector<uint8_t> _history, _stream;
    size_t _bits = 0, _bucket = 0, _fold_set = 0;
    uint32_t _c0 = 1;
    int _bit_count = 0, _expected = -1;
    uint8_t _last = 0;
    uint64_t _word = 0, _previous_word = 0;
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
// level (2). Then
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

// Analyzes the data, then compresses it the way that gives the smaller message: with or without the capitals
// transform for text, with the full or the bitmap model for folded data, unless the options impose one. The
// two ways are compared on the first sample bytes, then the whole data is compressed once the better way;
// with sample 0, or data that fits in the sample, both are compressed whole and the smaller message kept.
choice compress_best(const gmprime& file, const std::string& mode, const std::string& capitals_mode, const int forced_table_bits,
                     const size_t sample, const int level) {
    const auto data = file.bytes();
    const uint32_t check = crc32(data);
    choice best;
    best.cfg.level = level;

    // The first 128 KB find the same folds as the whole data, much sooner
    const auto folds = file.analyze_folds(64, 8192, 131072);
    best.cfg.fold = mode == "fold" || mode == "fold-full" || mode == "bitmap" || (mode == "auto" && folds.strong());
    if (best.cfg.fold) {
        best.cfg.width = folds.width;
        best.cfg.layer = folds.layer < 256 ? folds.layer : 0;
    }
    best.cfg.bitmap = mode == "bitmap";

    const auto lowered = capitals_encode(data);
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
        auto other = compare_capitals ? compress(first(*lowered), best.cfg, true, forced_table_bits, check)
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
    std::fprintf(stderr, "usage: gmz [-1|-2] [-m auto|text|fold|fold-full|bitmap] [-C auto|on|off] [-S bytes] [-T bits] c|d INPUT OUTPUT\n"
                         "       gmz [options] t INPUT\n");
    return 2;
}

} // namespace

int main(int argc, char** argv) {
    std::string mode = "auto", capitals_mode = "auto";
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
    const auto best = compress_best(input, mode, capitals_mode, forced_table_bits, sample, level);
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
    if (best.cfg.fold) {
        std::printf(" (width %zu bits, layer %zu%s)", best.cfg.width, best.cfg.layer, best.cfg.bitmap ? ", bitmap" : "");
    }
    std::printf("\n");
    return same ? 0 : 1;
}
