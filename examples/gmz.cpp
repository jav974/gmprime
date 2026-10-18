// gmz: a context-mixing compressor built on gmprime and its coding toolkit.
//
//   gmz [options] c INPUT OUTPUT    compress
//   gmz [options] d INPUT OUTPUT    decompress
//   gmz [options] t INPUT           compress, decompress, check and report, without writing anything
//
//   -m MODE   auto (default), text, fold (full and bitmap variants tried), fold-full, bitmap
//   -C MODE   capitals transform for text: auto (tried both ways, default), on, off
//   -T BITS   context table size in bits of slot count, instead of one sized for the data (the decompressor
//             needs the same value)
//
// The file is first analyzed with gmprime::analyze_folds. Data that folds (bitmap rows, fixed-size records)
// is coded with bit and byte templates across the fold, in 2D and 3D; anything else with the text model:
// context orders, words, line structure, last occurrences, indirect modeling, four mixers and chained maps.
// The compressed message is a gmprime: the arithmetic coder's output, saved and loaded like any number.
#include <chrono>
#include <cstdio>
#include <memory>
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
    int table_bits = 0;     // context tables of 2^table_bits slots
};

class model {
public:
    explicit model(const config& cfg): _cfg(cfg), _final(1, 1, 0.002), _map1(256 * 256, 0.02), _map2(size_t{1} << 16, 0.02), _match(6) {
        // Byte contexts, left out by the bitmap variant: orders with their counter limits (stable statistics for
        // low orders, adaptive ones above), words, line structure, last occurrences
        const bool bytes = !cfg.bitmap;
        if (bytes) {
            _orders = {0, 1, 2, 3, 4, 5, 6, 8};
            _limits = {255, 255, 127, 60, 30, 30, 30, 30};
        }
        _byte_contexts = _orders.size() + (bytes ? 2 + 2 + 1 + 1 + 1 + 1 : 0);
        _text = bytes;

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
        }
        _contexts = _byte_contexts + _bit_templates.size() + (cfg.fold ? 1 : 0) + _byte_templates;

        for (size_t t = 0; t < _contexts; ++t) {
            _tables.emplace_back(cfg.table_bits);
            _maps.emplace_back(1023);
            _limits.resize(std::max(_limits.size(), t + 1), 30);
        }
        _slots.resize(_contexts);
        _hashes.resize(_contexts);
        if (_text) {
            _last_word = std::make_unique<checked_table<last_byte_slot>>(cfg.table_bits - 2);
            _last_three = std::make_unique<checked_table<last_byte_slot>>(cfg.table_bits - 2);
        }

        // Mixers with weights selected by: the partial byte, the highest order seen with the match state, the
        // previous byte and the bit position, the previous byte's class, the fold neighbourhood
        const size_t inputs = _contexts * 3 + 2 + 1;
        _mixers.emplace_back(inputs, 256, 0.003);
        _mixers.emplace_back(inputs, (_orders.size() + 1) * 4, 0.003);
        _mixers.emplace_back(inputs, 256 * 8, 0.003);
        _mixers.emplace_back(inputs, 8 * 256, 0.003);
        if (cfg.fold) {
            _mixers.emplace_back(inputs, 256, 0.003);
        }
        _final = logistic_mixer(_mixers.size() + 1, 1, 0.002, 1.0 / static_cast<double>(_mixers.size()));
        _x.reserve(inputs);
        _stream.push_back(0);
        new_byte();
    }

    // Probability that the next bit is 1
    double p() {
        // Mixers learn faster at the start
        const double rate = 0.003 * (1 + 6.0 * 65536 / (65536 + static_cast<double>(_bits)));
        for (auto& mixer : _mixers) {
            mixer.set_rate(rate);
        }
        _x.clear();

        if (_cfg.fold) {
            fold_contexts();
        }

        // Each context: its probability, its history's probability across contexts, and a deterministic input
        // when only one bit value was ever seen
        size_t highest_seen = 0;
        for (size_t t = 0; t < _contexts; ++t) {
            const uint64_t key = t < _byte_contexts ? _hashes[t] ^ (static_cast<uint64_t>(_c0) * 0x9E3779B97F4A7C15ull) : _hashes[t];
            auto& slot = _tables[t].at(key);
            _slots[t] = &slot;
            _x.push_back(slot.seen() ? stretch(clamp_probability(slot.p)) : 0.0);
            _x.push_back(stretch(clamp_probability(_maps[t].p(slot.history))));
            const int zeros = slot.history.zeros(), ones = slot.history.ones();
            _x.push_back(zeros == 0 && ones > 0 ? ones / 4.0 : ones == 0 && zeros > 0 ? -zeros / 4.0 : 0.0);
            if (slot.seen() && t < _orders.size()) {
                highest_seen = t + 1;
            }
        }

        _expected = _match.expected_bit(_history, _c0, _bit_count);
        const uint32_t length = _match.length();
        _bucket = _expected < 0 ? 0 : std::min(length, 31u) + 1;
        const double sign = _expected == 1 ? 1.0 : -1.0;
        _x.push_back(_expected < 0 ? 0.0 : sign * stretch(clamp_probability(_match_counters[_bucket].p)));
        _x.push_back(_expected < 0 ? 0.0 : sign * std::min(length, 32u) / 16.0);
        const size_t match_state = _expected < 0 ? 0 : length < 16 ? 1 : length < 32 ? 2 : 3;
        _x.push_back(0.3);

        const size_t sets[] = {
            _c0,
            highest_seen * 4 + match_state,
            static_cast<size_t>(_last) * 8 + static_cast<size_t>(_bit_count),
            static_cast<size_t>(classify(_last)) * 256 + _c0,
            _fold_set,
        };
        for (size_t i = 0; i < _mixers.size(); ++i) {
            _mixers[i].add(_x);
            _final.add(stretch(_mixers[i].mix(sets[i])));
        }
        _final.add(0.3);
        _mixed = _final.mix(0);

        // Chained maps: order 1, then order 2 refining the first
        const double p1 = _map1.refine(_mixed, static_cast<size_t>(_last) * 256 + _c0);
        const double p2 = _map2.refine(p1, (hash_mix(context_hash(_history, _history.size(), 2)) ^ _c0) & 0xFFFF);
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
        _map1.update(bit);
        _map2.update(bit);

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
            _hashes[t++] = bit_template(_stream, pos, _bit_templates[k]) * 0x9E3779B97F4A7C15ull + 977 * (k + 1);
        }
        const auto runs = runs_at(_stream, pos, l);
        const uint64_t runs_key = (std::min<size_t>(runs.left, 63) << 6 | std::min<size_t>(runs.above, 63)) << 3 | bit_template(_stream, pos, {1, l, l - 1});
        _hashes[t++] = runs_key * 0x9E3779B97F4A7C15ull + 4444;

        if (_byte_templates > 0) {
            // The byte above (and two above, and one layer behind), aligned on the current byte's first bit
            const auto aligned = [&](const size_t distance) {
                std::vector<size_t> distances;
                for (size_t k = 0; k < 8; ++k) {
                    distances.push_back(distance + static_cast<size_t>(_bit_count) - k);
                }
                return bit_template(_stream, pos, distances);
            };
            const uint64_t above = aligned(l);
            const uint64_t above2 = aligned(2 * l);
            _hashes[t++] = (above << 9 | _c0) * 0x9E3779B97F4A7C15ull + 1111;
            _hashes[t++] = ((above2 << 8 | above) << 9 | _c0) * 0x9E3779B97F4A7C15ull + 2222;
            if (_byte_templates > 2) {
                const uint64_t behind = aligned(l * _cfg.layer);
                _hashes[t++] = ((behind << 8 | above) << 9 | _c0) * 0x9E3779B97F4A7C15ull + 3333;
            }
        }
        _fold_set = bit_template(_stream, pos, {1, 2, l - 1, l, l + 1, 2 * l, l - 2, l + 2});
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
        const uint64_t above = _lines.above(_history);
        _hashes[t++] = hash_mix(above * 257 + last) + 707;
        _hashes[t++] = hash_mix(above * 64 + std::min<size_t>(column, 63)) + 808;
        _hashes[t++] = hash_mix(std::min<size_t>(column, 80) * 1031 + std::min<size_t>(previous, 120) * 3 + (last == ' ')) + 909;
        const long to_end = static_cast<long>(column) - static_cast<long>(previous);
        _hashes[t++] = hash_mix(static_cast<uint64_t>(std::max(-40L, std::min(40L, to_end)) + 40) * 8
                                + static_cast<uint64_t>(classify(static_cast<uint8_t>(last)))) + 1010;

        // What followed the last occurrence of the current word, and of the last 3 bytes
        _word_slot = &_last_word->at(_word + 0x77);
        _hashes[t++] = hash_mix((_word_slot->count > 0 ? _word_slot->byte : 256u) * 257 + last) + 1212;
        _three_slot = &_last_three->at(context_hash(_history, pos, 3) + 0x99);
        _hashes[t++] = hash_mix((_three_slot->count > 0 ? _three_slot->byte : 256u) * 257 + last) + 1313;
    }

    const config& _cfg;
    std::vector<size_t> _orders;
    std::vector<int> _limits;
    size_t _byte_contexts = 0, _contexts = 0, _byte_templates = 0;
    bool _text = false;
    std::vector<std::vector<size_t>> _bit_templates;
    std::vector<checked_table<context_slot>> _tables;
    std::vector<history_map> _maps;
    std::vector<context_slot*> _slots;
    std::vector<uint64_t> _hashes;
    std::unique_ptr<checked_table<last_byte_slot>> _last_word, _last_three;
    last_byte_slot *_word_slot = nullptr, *_three_slot = nullptr;
    std::vector<logistic_mixer> _mixers;
    logistic_mixer _final;
    std::vector<double> _x;
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
    double _mixed = 0.5;
};

void encode_number(arithmetic_encoder& encoder, const size_t value, const int bits) {
    for (int i = bits - 1; i >= 0; --i) {
        encoder.encode(static_cast<int>((value >> i) & 1), 0.5);
    }
}

size_t decode_number(arithmetic_decoder& decoder, const int bits) {
    size_t value = 0;
    for (int i = 0; i < bits; ++i) {
        value = (value << 1) | static_cast<size_t>(decoder.decode(0.5));
    }
    return value;
}

// Header, at p = 1/2: payload size (32 bits), fold (1), bitmap (1), capitals (1), fold width (16), layer (8)
gmprime compress(const std::vector<uint8_t>& payload, config cfg, const bool capitals, const int forced_table_bits) {
    arithmetic_encoder encoder;
    encode_number(encoder, payload.size(), 32);
    encode_number(encoder, cfg.fold ? 1 : 0, 1);
    encode_number(encoder, cfg.bitmap ? 1 : 0, 1);
    encode_number(encoder, capitals ? 1 : 0, 1);
    encode_number(encoder, cfg.width, 16);
    encode_number(encoder, cfg.layer, 8);

    cfg.table_bits = forced_table_bits > 0 ? forced_table_bits : table_bits(payload.size());
    model m(cfg);
    for (const auto byte : payload) {
        for (int i = 7; i >= 0; --i) {
            const int bit = (byte >> i) & 1;
            encoder.encode(bit, m.p());
            m.update(bit);
        }
    }
    return encoder.finish();
}

std::vector<uint8_t> decompress(const gmprime& code, const int forced_table_bits) {
    arithmetic_decoder decoder(code);
    config cfg;
    const size_t size = decode_number(decoder, 32);
    cfg.fold = decode_number(decoder, 1) == 1;
    cfg.bitmap = decode_number(decoder, 1) == 1;
    const bool capitals = decode_number(decoder, 1) == 1;
    cfg.width = decode_number(decoder, 16);
    cfg.layer = decode_number(decoder, 8);

    cfg.table_bits = forced_table_bits > 0 ? forced_table_bits : table_bits(size);
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
    return capitals ? capitals_decode(payload) : payload;
}

struct choice {
    gmprime code;
    config cfg;
    bool capitals = false;
};

// Analyzes the data, then compresses it the way that gives the smallest message among those tried
choice compress_best(const gmprime& file, const std::string& mode, const std::string& capitals_mode, const int forced_table_bits) {
    const auto data = file.bytes();
    choice best;

    const auto folds = file.analyze_folds();
    best.cfg.fold = mode == "fold" || mode == "fold-full" || mode == "bitmap" || (mode == "auto" && folds.strong());
    if (best.cfg.fold) {
        best.cfg.width = folds.width;
        best.cfg.layer = folds.layer < 256 ? folds.layer : 0;
    }

    // Capitals transform: for text only, tried both ways unless told otherwise
    const auto lowered = capitals_encode(data);
    const bool try_capitals = !best.cfg.fold && lowered.has_value() && capitals_mode != "off";
    if (capitals_mode != "on" || !try_capitals) {
        best.code = compress(data, best.cfg, false, forced_table_bits);
    }

    // Folded data: the bitmap variant too, the smaller kept
    if (best.cfg.fold && mode != "fold-full") {
        config variant = best.cfg;
        variant.bitmap = true;
        auto code = compress(data, variant, false, forced_table_bits);
        if (mode == "bitmap" || code.bytes().size() < best.code.bytes().size()) {
            best.code = std::move(code);
            best.cfg = variant;
        }
    }

    if (try_capitals) {
        auto code = compress(*lowered, best.cfg, true, forced_table_bits);
        if (capitals_mode == "on" || code.bytes().size() < best.code.bytes().size()) {
            best.code = std::move(code);
            best.capitals = true;
        }
    }
    return best;
}

int usage() {
    std::fprintf(stderr, "usage: gmz [-m auto|text|fold|fold-full|bitmap] [-C auto|on|off] [-T bits] c|d INPUT OUTPUT\n"
                         "       gmz [options] t INPUT\n");
    return 2;
}

} // namespace

int main(int argc, char** argv) {
    std::string mode = "auto", capitals_mode = "auto";
    int forced_table_bits = 0;
    int arg = 1;
    for (; arg + 1 < argc && argv[arg][0] == '-'; arg += 2) {
        const std::string flag = argv[arg];
        if (flag == "-m") mode = argv[arg + 1];
        else if (flag == "-C") capitals_mode = argv[arg + 1];
        else if (flag == "-T") forced_table_bits = std::stoi(argv[arg + 1]);
        else return usage();
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
        gmprime::from_bytes(decompress(input, forced_table_bits)).save_to_file(argv[arg + 2]);
        return 0;
    }

    const auto start = std::chrono::steady_clock::now();
    const auto best = compress_best(input, mode, capitals_mode, forced_table_bits);
    if (command == "c") {
        best.code.save_to_file(argv[arg + 2]);
        return 0;
    }

    // Test: decode what was coded and report
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const auto data = input.bytes();
    const bool same = decompress(best.code, forced_table_bits) == data;
    const size_t size = best.code.bytes().size();
    std::printf("%zu -> %zu bytes, %.3f bits/byte, %s, %.1f s, %s%s", data.size(), size, data.empty() ? 0.0 : 8.0 * size / data.size(),
                same ? "round trip ok" : "ROUND TRIP FAILED", seconds, best.cfg.fold ? "fold" : "text", best.capitals ? " + capitals" : "");
    if (best.cfg.fold) {
        std::printf(" (width %zu bits, layer %zu%s)", best.cfg.width, best.cfg.layer, best.cfg.bitmap ? ", bitmap" : "");
    }
    std::printf("\n");
    return same ? 0 : 1;
}
