![Build Status](https://github.com/jav974/gmprime/actions/workflows/cmake-multi-platform.yml/badge.svg)
![License](https://img.shields.io/github/license/jav974/gmprime)
![Coverage](https://img.shields.io/codecov/c/github/jav974/gmprime)

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
