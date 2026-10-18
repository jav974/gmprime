// Writes synthetic files with the structures gmz looks for, for its tests:
//   records.bin   records of 28 bytes: an identifier, decimal doubles of either byte order (a lattice), an IBM
//                 single precision value and a float, with a few values off their lattice
//   image_le.bin  a smooth image of 16-bit little-endian pixels, rows of 72 pixels
//   image_be.bin  the same in big-endian, rows of 60 pixels
// Each just large enough for gmz to find its structure, so that the tests stay short
//   gmz_samples DIRECTORY
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace {

uint32_t state = 12345;

uint32_t next() {
    state = state * 1664525u + 1013904223u;
    return state >> 8;
}

void put(std::vector<uint8_t>& out, const uint64_t bits, const int size, const bool big_endian) {
    for (int i = 0; i < size; ++i) {
        const int shift = 8 * (big_endian ? size - 1 - i : i);
        out.push_back(static_cast<uint8_t>(bits >> shift));
    }
}

uint64_t double_bits(const double v) {
    uint64_t b;
    std::memcpy(&b, &v, 8);
    return b;
}

uint32_t float_bits(const float v) {
    uint32_t b;
    std::memcpy(&b, &v, 4);
    return b;
}

// IBM single precision of k / 64 for 0 < k < 2^24
uint32_t ibm_bits(uint32_t k) {
    int exponent = 64 + 6;    // fraction = k / 2^24 * 16^6, times 16^(-6) / 64... built by normalizing below
    uint64_t fraction = uint64_t{k} << 18;    // k / 64 = k * 2^18 / 2^24
    exponent = 64;
    while (fraction >= (uint64_t{1} << 24)) {
        fraction >>= 4;
        ++exponent;
    }
    while (fraction != 0 && fraction < (uint64_t{1} << 20)) {
        fraction <<= 4;
        --exponent;
    }
    return static_cast<uint32_t>(exponent) << 24 | static_cast<uint32_t>(fraction);
}

void save(const std::string& path, const std::vector<uint8_t>& bytes) {
    std::ofstream(path, std::ios::binary).write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
}

std::vector<uint8_t> image(const size_t width, const size_t height, const bool big_endian) {
    std::vector<uint8_t> out;
    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            const auto value = static_cast<uint64_t>(1000 + 3 * x + 5 * y + (x * y) % 7 + next() % 4);
            put(out, value, 2, big_endian);
        }
    }
    return out;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        return 2;
    }
    const std::string dir = std::string(argv[1]) + "/";

    std::vector<uint8_t> records;
    int64_t ra = 100000, dec = -5000;
    for (uint32_t r = 0; r < 300; ++r) {
        ra += next() % 50;
        dec += static_cast<int64_t>(next() % 21) - 10;
        put(records, r, 4, false);
        put(records, double_bits(r % 97 == 5 ? 1.0 / 3 : static_cast<double>(ra) * 0.001), 8, false);
        put(records, double_bits(static_cast<double>(dec) * 0.01), 8, true);
        put(records, ibm_bits(1 + next() % 5000), 4, true);
        put(records, float_bits(static_cast<float>(next() % 300) * 0.25f), 4, false);
    }
    save(dir + "records.bin", records);
    save(dir + "image_le.bin", image(72, 64, false));
    save(dir + "image_be.bin", image(60, 80, true));
    return 0;
}
