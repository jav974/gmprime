![Build Status](https://github.com/jav974/gmprime/actions/workflows/cmake-multi-platform.yml/badge.svg)
[![codecov](https://codecov.io/gh/jav974/gmprime/graph/badge.svg?token=QM7E5GTC4S)](https://codecov.io/gh/jav974/gmprime)
![License](https://img.shields.io/github/license/jav974/gmprime)

# GMPrime

GMPrime is a C++ library that provides a robust wrapper for the GNU Multiple Precision Arithmetic Library (GMP). It allows for deep integration with C++, making it easier to perform high-precision arithmetic operations.

## Features

- **Arithmetic Operations**: Supports all builtin arithmetic operations including addition, subtraction, multiplication, division, and modulo.
- **Bitwise Operations**: Includes bitwise operations such as AND, OR, XOR, and shifts.
- **Comparison Operations**: Enables comparison operations like equality, inequality, and relational comparisons.
- **File Operations**: Allows loading from and saving to files.
- **Utilities**: Includes utility functions for prime checking, factorial calculation, and more.
- **Byte/Bit access**: Enables direct read/write access to bytes and bits within the gmprime object through the `operator[]`.
- **Custom Literals**: Provides custom literals for intuitive byte/bit access.

## Quick start

### Prerequisites

- **CMake**: Ensure CMake is installed on your system.
- **C++ Compiler**: Ensure you have a C++ compiler that supports C++17 or later.
- **GMP Library**: Make sure the GMP library is installed on your system.

### CMake Integration

To use GMPrime in your CMake project, add the following to your `CMakeLists.txt`:

```cmake
include(FetchContent)
FetchContent_Declare(
        gmprime
        GIT_REPOSITORY https://github.com/jav974/gmprime.git
        GIT_TAG main
)
FetchContent_MakeAvailable(gmprime)
target_link_libraries(${PROJECT_NAME} gmprime)
```

### Build from source (optional)

1. Clone the repository:
    ```sh
    git clone https://github.com/jav974/gmprime.git
    cd gmprime
    ```

2. Create a build directory and run CMake:
    ```sh
    mkdir build
    cd build
    cmake ..
    make
    ```
   *This will create the gmprime static library, as well as an executable for tests.*

3. Run tests:
    ```sh
   ctest --output-on-failure
    ```

### Linux test matrix (Docker)

Builds and tests with GCC and Clang, on x86_64 and ARM64, in Release then in Debug with AddressSanitizer and UndefinedBehaviorSanitizer:

```sh
docker/test-all.sh                  # every combination in parallel, with a summary
docker compose run --rm gcc-amd64   # a single one (see compose.yaml for the list)
```

Sources are mounted read-only and build trees are kept in Docker volumes, so later runs only rebuild what changed. The non-native architecture runs emulated. `docker compose down -v` removes the build volumes.

Before pushing, `docker/coverage.sh` measures coverage like the CI and lists the lines added since `main` (or the revision given as argument), committed or not, that no test executes: what Codecov reports as patch coverage.

## Basic usage

```cpp
#include <iostream>
#include <gmprime.hpp>

using namespace GMPrime;

int main() {
    gmprime a = 123;
    gmprime b("2A", 16);

    auto sum = a + b;

    std::cout << "Sum: " << sum << std::endl;
    return 0;
}
```

### Bit windows

`bits(pos, width)` accesses `width` bits starting at bit `pos`, at any position and of any width. A window holds bits of the magnitude: a sign has no meaning in it, so only unsigned integers and `gmprime` convert to and from it.

```cpp
gmprime n("0x0123456789abcdef", 0);

uint8_t byte = n.bits(4, 8);     // 0xde; std::overflow_error if the window is wider than the type
gmprime wide = n.bits(8, 40);    // 0x456789abcd, any width
n.bits(0, 4) = 0u;               // keeps the low 4 bits of the value, like a C bit field

std::cout << gmprime::hex << n;  // 123456789abcde0
```

### Bit analysis

Counting and scanning work on the magnitude like every bit access, a limb at a time:

```cpp
gmprime text;
text.load_from_file("data.txt");

size_t ones = text.popcount();                    // number of bits set
size_t diff = text.hamming_distance(other);       // number of bits that differ
size_t first = text.next_set_bit();               // gmprime::npos if there is none
size_t hole = text.next_clear_bit(first);         // always found: bits beyond the magnitude are clear

for (const size_t pos : (text & mask).set_bits()) {
    // positions of the bits set, in increasing order: whole zero limbs are skipped
}
```

Symbol statistics count the values of `width`-bit windows, aligned or sliding, over the whole number in whole bytes, like a file (a loaded file's high zero bytes included):

```cpp
std::vector<size_t> bytes = text.histogram(8);       // bytes[c]: occurrences of the byte c
std::vector<size_t> sliding = text.histogram(8, 1);  // 8-bit windows starting at every bit
double bits_per_byte = gmprime::entropy(bytes);      // Shannon entropy, in bits per symbol
```

The entropy of a histogram is a lower bound for any code encoding its symbols independently of each other. Compressors that model context or repetitions can go below it.

### Context and repetitions

The two kinds of structure compressors exploit can be measured on their own: the context of each symbol (which byte tends to follow which), and repeats of whole strings.

```cpp
auto order2 = text.conditional_entropy(2);         // entropy of a byte given the 2 before it
double naive = order2.naive;                       // bits per byte, optimistic when contexts are rare
double corrected = order2.corrected;               // with the Miller-Madow correction
double bits = text.adaptive_code_length(2);        // achievable size of an adaptive order-2 model

auto lpf = text.longest_previous_factors();        // lpf[i]: longest string at byte i seen before
auto parse = text.repetitions(4);                  // greedy LZ77 parse with repeats of 4 bytes or more
size_t covered = parse.matched_bytes;              // also parse.matches, parse.literals
```

Entropy estimates keep decreasing with the order, as contexts get too rare to be measured reliably: `adaptive_code_length` is the honest figure, the bits a context model actually takes to code the data while learning it, comparable to the output of a real compressor. Repetitions are what dictionary compressors (LZ77, xz) reference, and no context model sees them.

### Folds

Fold a bit field at width L, like fabric, and bits L apart stack up. `autocorrelation` measures, for every distance, how many bits equal the one that far ahead, with one shift, XOR and popcount of the whole number per distance:

```cpp
auto self = data.autocorrelation(4096);
// Peaks of self[d] / (self[0] - d) are fold widths: 1728 bits for the rows of a fax image,
// 104 bits for 13-byte records, about a line length for formatted text
```

`analyze_folds` does the search for you, in stream order (each byte from its highest bit, the order of a bitmap's pixels, which `stream_order()` gives to any number), and scores each width against independent bits:

```cpp
auto folds = data.analyze_folds();  // widths from 64 to 8192 bits; analyze_folds(64, 8192, 131072)
                                    // looks at the first 128 KB only, much sooner on large data
if (folds.strong()) {               // the best width stands out: bitmap rows, fixed-size records
    folds.width;                    // bits per row
    folds.layer;                    // rows per layer, for a 3D view (0: none)
}
```

To build a compressor from these measurements, see the [coding toolkit](#coding-toolkit).

### Enumerative coding

A bit field of `n` bits with `k` of them set is one of `C(n, k)` arrangements. Its rank among them is a number below `C(n, k)`: `log2 C(n, k)` bits, exactly the information left once `n` and `k` are known.

```cpp
gmprime mask = ...;                                       // 10,000 bits, 100 of them set
gmprime rank = mask.combination_rank(10'000);             // 804 bits, where the positions take 1,400
gmprime same = gmprime::from_combination_rank(rank, 10'000, 100);
gmprime count = gmprime::binomial(10'000, 100);           // the number of arrangements
```

The length is given because leading zeros are part of the field. Ranking takes time quadratic in the length: for long fields, `rank_blocks` in the [coding toolkit](#enumerative-coding-1) ranks blocks and packs their ranks into one number. The rank is the shortest code for bits whose only regularity is their count, such as a random sparse mask. Data with more structure, like text or images, is better served by context models.

### Integer behavior

A `gmprime` is used like an `int`: arithmetic, bitwise and comparison operators (with a primitive on either side), `if (a)`, streams with every formatting flag (`std::showbase`, `std::setw`...) and `>>`, `abs`, `gcd`, `lcm`, `to_string`, `std::hash` and `std::numeric_limits`. Comparisons with a floating point value are exact: `gmprime(5) < 5.5`.

Where an `int` has undefined behavior or overflows, a `gmprime` throws: `std::domain_error` on division by zero, a negative shift count or the square root of a negative number, `std::length_error` on a left shift beyond what GMP can hold.

The differences are deliberate:

- Only the prefix forms of `++` and `--` exist, to avoid copying a large number for nothing.
- Bit access (`[]`, `get_bit`, `bits`, `as_span`) reads the magnitude as GMP stores it, while bitwise operators (`&`, `|`, `^`, `~`, `>>`) follow two's complement like an `int`. For -2, `get_bit(2)` is 0 but `((a >> 2) & 1)` is 1.
- A negative number prints as a sign and its magnitude in every base (`-0xff`), where an `int` prints its two's complement in octal and hexadecimal.
- There is no implicit conversion to a primitive, which could truncate: use `to_primitive<T>()`.
- Mixed arithmetic with a floating point value converts it to `gmprime` first, truncating it: `a + 0.5` is `a`. Use `a.to_primitive<double>() + 0.5` for a floating point result.

## Coding toolkit

`#include <gmprime_coding.hpp>` adds the tools to build a compressor on gmprime: models that predict each bit of the data from what precedes it, and an arithmetic coder that turns those predictions into a compressed message, itself a `gmprime` to save, load or analyze like any other.

### Principles

- **Bit by bit.** Each byte is coded from its highest bit. While coding a byte, `c0` holds the bits coded so far behind a leading 1 (1 at the start of the byte, up to 255), and `bit_count` how many there are. Contexts usually combine the previous bytes with `c0`.
- **A probability per bit.** A model gives the probability `p` that the next bit is 1, on 16 bits (0 to 65535 for 0 to 1). The coder spends about `-log2(p)` bits on a 1 and `-log2(1 - p)` on a 0, so good predictions mean a small message. `clamp_probability` keeps `p` within what the coder handles.
- **Causal models, identical on both sides.** For every bit: ask the model for `p`, code (or decode) the bit, then let the model learn it. Models only use what precedes the bit, so the decoder, running the same code on the bits it has decoded, makes the same predictions.
- **Integers only.** Models compute in fixed point, without floating point: the same model gives the same predictions with any compiler on any platform, so a message coded on one decodes on all the others. Floating point does not guarantee it (fused multiply-add, for example, changes the last bits of a result, and a single different prediction corrupts the rest of the message).

### Probabilities

| Tool | Role |
|---|---|
| `stretch(p)`, `squash(x)` | the logistic domain, where predictions are mixed: `stretch` is `ln(p / (1 - p))` in 1/256 units (-3071 to 3071), `squash` its inverse, interpolated from a table and clamped beyond ±12 |
| `clamp_probability(p)`, `probability_12(p)` | keep `p` within the coder's 12 bits |
| `bit_counter` | adaptive probability: `update(bit, limit)` moves it toward the bit at rate `1 / (n + 1.5)`, fast at first, then more slowly as `n` grows up to `limit` (high for stable data, low to follow changes); `p16()` gives it on 16 bits (it keeps 22, so that slow rates still move it); `seen()` tells whether it learned anything |
| `bit_history` | the counts of 0s and 1s seen in a context in one byte (`zeros()`, `ones()`, up to 15), recent bits weighing more |
| `history_map` | indirect modeling: `p(history)` is the probability learned, across every context, for contexts with that history; `update(history, bit)` after coding. Rare contexts know little alone, their history generalizes |

### Contexts

| Tool | Role |
|---|---|
| `context_hash(history, pos, order)` | hash of the `order` bytes before `pos`, distinct per order |
| `hashed_table<T>(bits)` | `2^bits` slots of type `T`; `at(hash)` returns the slot of a context, colliding contexts sharing it; `hashed_counters` holds `bit_counter`s |
| `checked_table<T>(bits)` | the same with an 8-bit checksum in each slot: a context landing on another's slot resets it rather than sharing its statistics (`T` needs a `check` field). Created from zeroed memory, so only touched pages are used. `prefetch(hash)` brings a slot into cache ahead of time, e.g. both possible slots of the next bit; `at(hash, offset)` addresses 16 adjacent slots per hash, e.g. the nodes of a nibble |
| `context_slot` | 8 bytes per context: a probability (as `bit_counter`), its `bit_history` and the checksum; `update(bit, limit)` updates both |
| `last_byte_slot` | the byte that followed the context last time and how many times in a row (`update(byte)`); `expected_bit(c0, bit_count)` predicts its next bit, or -1 if the byte coded so far disagrees |

A table takes `2^bits` times the slot size: 32 MB for `checked_table<context_slot>(22)`.

### Mixing

| Tool | Role |
|---|---|
| `logistic_mixer(inputs, sets, rate)` | combines predictions: `add(x)` each input (usually `stretch(p)`, 0 for "no idea"), `mix(set)` returns `p`, `update(bit)` learns the weights. Each set has its own weights, chosen by a small context (e.g. `c0`), so the mixer learns whom to trust when. Rates are scaled by 65536 (197 for 0.003); `set_rate` changes it, e.g. decreasing over time. `mix(inputs, set)` reads a `std::vector<int>` kept by the caller instead: several mixers can share the same inputs without copying them |
| `probability_map(contexts, rate)` | refines a probability given a context: `refine(p, context)`, then `update(bit)`; rate scaled by 65536 (1311 for 0.02). Averaged with its input, it corrects what the mixer systematically gets wrong |
| `match_model(min_length)` | finds the last earlier occurrence of the last `min_length` bytes and predicts that what followed comes again: `byte_added(history)` after each byte, `expected_bit(history, c0, bit_count)`, `length()` of the match, `miss()` when a prediction fails |

### Text and bit geometry

| Tool | Role |
|---|---|
| `classify(byte)` | the `char_class` of a byte (newline, space, lowercase, uppercase, digit, punctuation, other), a good mixer set or context for text |
| `line_tracker` | the 2D view of a text: `column(history)`, `above(history)` (byte at the same column of the previous line), `previous_length()`; `byte_added(history)` after each byte |
| `capitals_encode(text)`, `capitals_decode(encoded)` | reversible lowercasing: a capitalized word is lowercased behind `0x01`, an all-capitals word behind `0x02`, so that "The" and "the" share statistics; `capitals_encode` returns nothing if the text holds those bytes |
| `stream_bit(bytes, p)` | bit `p` of bytes read as a stream, each byte from its highest bit, the order of a bitmap's pixels |
| `bit_template(bytes, pos, distances)` | the stream bits at `pos - d` for each distance `d`, packed into a context. With the fold width `L` found by `gmprime::analyze_folds`, distances such as 1, 2, `L - 1`, `L`, `L + 1` form a 2D neighbourhood, and `L * H` reaches the layer behind in 3D |
| `stream_bits(bytes, first, count)` | up to 57 consecutive stream bits in one word read, the first one highest; bits before the start (`first` may be negative) and beyond the end are 0 |
| `compiled_template(distances)` | `bit_template` prepared once: `at(bytes, pos)` gives the same context, reading each run of consecutive distances (`L + 1, L, L - 1`...) as one word. For templates used at every bit |
| `runs_at(bytes, pos, width)` | pixel runs of a bitmap folded at `width` bits: `left`, the run of equal bits ending before `pos`, and `above`, how far the bit above stays the same to the right |

### Enumerative coding

| Tool | Role |
|---|---|
| `rank_blocks(bits, length, block_bits)`, `unrank_blocks(ones, rank, length, block_bits)` | a bit field ranked by blocks: the bits set in each block (`ones`, for a model to code) and the blocks' combination ranks packed into one number below `radix`, the product of their binomials. Linear time, and the counts per block follow the density where it changes |
| `mixed_radix_pack(digits, radices)`, `mixed_radix_unpack(value, radices)` | numbers in different bases, each below its own radix, as one number: no fraction of a bit lost per number. Over a product tree, so packing many is fast |
| `multiset_rank(bytes)`, `from_multiset_rank(rank, counts)`, `multiset_arrangements(counts)` | the same for bytes: once the count of each byte value is known, the rank of their order among the `n! / (product of count!)` arrangements |

### Coding and measuring

| Tool | Role |
|---|---|
| `arithmetic_encoder` | `encode(bit, p)` for each bit, then `finish()` returns the message as a `gmprime` |
| `arithmetic_decoder(code)` | `decode(p)` returns each bit, given the same `p` as the encoder. The decoder must know when to stop: store the length in the message |
| `cost_profile` | `add(bit, p)` alongside the coder: `per_byte()` gives the bits each byte cost, `total()` the size, to see where a model pays |
| `crc32(bytes)` | the CRC-32 of bytes (as zip and png): stored with a message, it tells whether decoding gave back the original data |
| `table_bits(data_bytes)` | a context table size for the data: about four slots per bit, from 16 to 22 bits, so small inputs take little memory; compressor and decompressor derive it from the stored length |

### A complete compressor

A model mixing orders 1 to 4 (with their bit histories) and a match model, shared by the compressor and the decompressor:

```cpp
#include <gmprime_coding.hpp>

using namespace GMPrime;
using namespace GMPrime::coding;

class model {
public:
    model(): _mixer(2 * 4 + 2 + 1, 256, 328), _map(256 * 256, 1311), _match(6) {
        for (int i = 0; i < 4; ++i) {
            _tables.emplace_back(22);
        }
    }

    // Probability that the next bit is 1, on 16 bits
    int p() {
        for (size_t i = 0; i < _tables.size(); ++i) {
            _slots[i] = &_tables[i].at(context_hash(_history, _history.size(), i + 1) * 256 + _c0);
            _mixer.add(_slots[i]->seen() ? stretch(clamp_probability(_slots[i]->p16())) : 0);
            _mixer.add(stretch(clamp_probability(_maps[i].p(_slots[i]->history))));
        }
        _expected = _match.expected_bit(_history, _c0, _bit_count);
        const int sign = _expected == 1 ? 1 : -1;
        _mixer.add(_expected < 0 ? 0 : sign * static_cast<int>(std::min(_match.length(), 32u)) * 64);
        _mixer.add(_expected < 0 ? 0 : sign * 256);
        _mixer.add(77);  // bias

        const int mixed = _mixer.mix(_c0);
        const uint8_t previous = _history.empty() ? 0 : _history.back();
        return clamp_probability((mixed + 3 * _map.refine(mixed, previous * 256u + _c0)) / 4);
    }

    // Learns the bit, after it was coded or decoded
    void update(const int bit) {
        _mixer.update(bit);
        _map.update(bit);
        for (size_t i = 0; i < _tables.size(); ++i) {
            _maps[i].update(_slots[i]->history, bit);
            _slots[i]->update(bit, 60);
        }
        if (_expected >= 0 && bit != _expected) {
            _match.miss();
        }
        _c0 = (_c0 << 1) | static_cast<uint32_t>(bit);
        if (++_bit_count == 8) {
            _history.push_back(static_cast<uint8_t>(_c0));
            _match.byte_added(_history);
            _c0 = 1;
            _bit_count = 0;
        }
    }

private:
    std::vector<checked_table<context_slot>> _tables;
    history_map _maps[4];
    context_slot* _slots[4] = {};
    logistic_mixer _mixer;
    probability_map _map;
    match_model _match;
    std::vector<uint8_t> _history;
    uint32_t _c0 = 1;
    int _bit_count = 0, _expected = -1;
};

gmprime compress(const std::vector<uint8_t>& data) {
    arithmetic_encoder encoder;
    for (int i = 31; i >= 0; --i) {
        encoder.encode(static_cast<int>(data.size() >> i) & 1, 32768);  // the length, for the decoder
    }
    model m;
    for (const uint8_t byte : data) {
        for (int i = 7; i >= 0; --i) {
            const int bit = (byte >> i) & 1;
            encoder.encode(bit, m.p());
            m.update(bit);
        }
    }
    return encoder.finish();
}

std::vector<uint8_t> decompress(const gmprime& code) {
    arithmetic_decoder decoder(code);
    size_t size = 0;
    for (int i = 0; i < 32; ++i) {
        size = (size << 1) | static_cast<size_t>(decoder.decode(32768));
    }
    model m;
    std::vector<uint8_t> data;
    while (data.size() < size) {
        int byte = 0;
        for (int i = 0; i < 8; ++i) {
            const int bit = decoder.decode(m.p());
            m.update(bit);
            byte = (byte << 1) | bit;
        }
        data.push_back(static_cast<uint8_t>(byte));
    }
    return data;
}
```

Used with files:

```cpp
gmprime file;
file.load_from_file("alice29.txt");
compress(file.bytes()).save_to_file("alice29.gmc");

gmprime code;
code.load_from_file("alice29.gmc");
std::vector<uint8_t> data = decompress(code);
```

This model compresses `alice29.txt` (152,089 bytes, from the Canterbury corpus) to 38,871 bytes, where bzip2 -9 gives 43,202 and xz -9e 48,536. A research compressor adding the text tools, more orders and a second mixer reaches 1.95 bits per character on it, and folding at the width found by `autocorrelation` halves what xz achieves on a fax image.

### Going further

- **More contexts, more mixers.** Add orders (6, 8), a word context (hash of the current word's letters), `line_tracker` contexts for text; select mixer weights by `c0`, by the previous byte's `classify`, or by the highest order already seen, and mix several mixers with a final one.
- **Look at the data first.** `gmprime::analyze_folds` reveals fold widths (rows of an image, fixed-size records): for such data, contexts from `bit_template` across the fold beat any byte context.
- **Measure.** `cost_profile` shows which bytes cost the most, and `adaptive_code_length` or `conditional_entropy` how much structure a context order holds.

`examples/gmz.cpp` puts all of this together.

## Example: gmz

`gmz` is a compressor built on gmprime and its coding toolkit, in `examples/` (built with the project, `GMPRIME_BUILD_EXAMPLES` to turn it off):

```sh
gmz c data.txt data.gmz    # compress
gmz d data.gmz data.txt    # decompress
gmz t data.txt             # compress, decompress, check, report
gmz -1 c data.txt data.gmz # faster, a little larger (-2: a little more of both)
```

It analyzes the file first with `analyze_folds`. Data that folds is coded with bit and byte templates across the fold, in 2D and 3D, and pixel runs; anything else with the text model: orders 0 to 8, words, line structure, last occurrences, indirect modeling, four mixers selected by the partial byte, the highest order seen, the previous byte and its class, and chained probability maps. It tries the capitals transform on text, and a bitmap variant on folded data, and keeps the smaller message. The two ways are compared on the first 64 KB, then the data is compressed once. A CRC-32 of the original ends the message: decompressing anything else (a corrupted file, another table size than the compression's) fails instead of writing wrong data. Options: `-1` and `-2` choose a faster level (below), `-m` forces a mode, `-C` the capitals transform, `-S` the sample (0 compresses the data both ways, slower and sometimes a little smaller), `-T` a table size.

On the Canterbury corpus, sizes in bytes:

| File | Size | gzip -9 | bzip2 -9 | xz -9e | gmz | Mode found |
|---|---|---|---|---|---|---|
| alice29.txt | 152,089 | 54,194 | 43,202 | 48,536 | **35,567** | text + capitals |
| asyoulik.txt | 125,179 | 48,803 | 39,569 | 44,600 | **33,545** | text + capitals |
| lcet10.txt | 426,754 | 144,462 | 107,706 | 119,496 | **86,152** | text + capitals |
| plrabn12.txt | 481,861 | 194,357 | 145,577 | 165,464 | **123,387** | text + capitals |
| ptt5 (fax image) | 513,216 | 52,238 | 49,759 | 39,868 | **21,863** | fold of 1728 bits |
| kennedy.xls | 1,029,744 | 207,053 | 130,280 | 51,876 | **12,969** | fold of 104 bits, layers of 4 |
| cp.html | 24,603 | 7,960 | 7,624 | 7,660 | **6,177** | text + capitals |
| fields.c | 11,150 | 3,136 | 3,039 | 3,040 | **2,323** | text + capitals |
| grammar.lsp | 3,721 | 1,246 | 1,283 | 1,300 | **994** | text |
| sum | 38,240 | 12,854 | 12,909 | 9,508 | **9,413** | text |
| xargs.1 | 4,227 | 1,756 | 1,762 | 1,820 | **1,444** | text + capitals |

The price is speed, as for every context-mixing compressor. On Apple Silicon:

| File | gmz compress | gmz decompress | xz -9e compress | xz decompress |
|---|---|---|---|---|
| alice29.txt | 1.1 s | 0.5 s | 0.05 s | 0.01 s |
| lcet10.txt | 1.9 s | 1.2 s | 0.08 s | 0.01 s |
| ptt5 | 1.3 s | 0.8 s | 0.15 s | 0.01 s |
| kennedy.xls | 4.6 s | 4.1 s | 0.40 s | 0.02 s |

The faster levels trade a little size for speed (compress / decompress):

| Level | Canterbury total | alice29.txt | lcet10.txt | kennedy.xls |
|---|---|---|---|---|
| default | 333,834 bytes | 1.1 s / 0.5 s | 1.9 s / 1.2 s | 4.6 s / 4.1 s |
| `-1` | 337,489 bytes (+1.1%) | 0.4 s / 0.3 s | 0.9 s / 0.8 s | 3.2 s / 2.9 s |
| `-2` | 339,839 bytes (+1.8%) | 0.4 s / 0.2 s | 0.9 s / 0.7 s | 3.0 s / 2.6 s |

For text they keep fewer contexts and mixers, group the slots of each half byte, use smaller tables and skip the capitals comparison; folded data, whose time goes to the fold templates, only loses its text contexts. `-2` also drops the chained maps. The level is stored in the message.

Messages are portable: the same file gives the same message with GCC and Clang, on x86_64 and ARM64, Linux and macOS, and each decodes the others'.

Memory reaches about 800 MB for large inputs (tables are sized for the data, so small files take much less). Decompression runs the same model once, so it takes about as long as compressing does after the comparison. How these times were reached, 2.5 to 6 times faster than the first version, is told in [docs/RESEARCH.md](docs/RESEARCH.md).
