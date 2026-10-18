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

### Integer behavior

A `gmprime` is used like an `int`: arithmetic, bitwise and comparison operators (with a primitive on either side), `if (a)`, streams with every formatting flag (`std::showbase`, `std::setw`...) and `>>`, `abs`, `gcd`, `lcm`, `to_string`, `std::hash` and `std::numeric_limits`. Comparisons with a floating point value are exact: `gmprime(5) < 5.5`.

Where an `int` has undefined behavior or overflows, a `gmprime` throws: `std::domain_error` on division by zero, a negative shift count or the square root of a negative number, `std::length_error` on a left shift beyond what GMP can hold.

The differences are deliberate:

- Only the prefix forms of `++` and `--` exist, to avoid copying a large number for nothing.
- Bit access (`[]`, `get_bit`, `bits`, `as_span`) reads the magnitude as GMP stores it, while bitwise operators (`&`, `|`, `^`, `~`, `>>`) follow two's complement like an `int`. For -2, `get_bit(2)` is 0 but `((a >> 2) & 1)` is 1.
- A negative number prints as a sign and its magnitude in every base (`-0xff`), where an `int` prints its two's complement in octal and hexadecimal.
- There is no implicit conversion to a primitive, which could truncate: use `to_primitive<T>()`.
- Mixed arithmetic with a floating point value converts it to `gmprime` first, truncating it: `a + 0.5` is `a`. Use `a.to_primitive<double>() + 0.5` for a floating point result.
