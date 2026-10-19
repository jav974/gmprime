// Writes synthetic files with the structures gmz looks for, for its tests:
//   records.bin   records of 28 bytes: an identifier, decimal doubles of either byte order (a lattice), an IBM
//                 single precision value and a float, with a few values off their lattice
//   image_le.bin  a smooth image of 16-bit little-endian pixels, rows of 72 pixels
//   image_be.bin  the same in big-endian, rows of 60 pixels
//   numbers.csv   a time series: a count of seconds, the same time as a local date (UTC-7), a price printed from a
//                 float32 with 6 decimals, a temperature, a counter, the time in milliseconds, a quantity and the
//                 price times it, a smooth series
//   oneline.json  records as JSON on a single line
//   numbers.json  records as JSON, keyed numbers of various precisions, a title repeating one of them
//   numbers.log   log lines with dates, addresses, counters, random identifiers and hexadecimal
//   crowded.txt   numbers, and every byte that could mark them, so that none can
//   polylines.bin shapes of varying length, as GIS files store them: a count, then points as decimal doubles
//   header_records.bin  a header, then records of 3 decimal doubles
// Each just large enough for gmz to find its structure, so that the tests stay short
//   gmz_samples DIRECTORY
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
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

std::vector<uint8_t> text(const std::string& s) {
    return {s.begin(), s.end()};
}

std::string time_series() {
    std::string out = "epoch,local,price,temperature,count,millis,quantity,total,square\n";
    int64_t seconds = 1128000000;
    float price = 127.49f;
    int temperature = 215, count = 1000;
    char line[160];
    for (int i = 0; i < 160; ++i) {
        seconds += 60 + next() % 3;
        const std::time_t local = static_cast<std::time_t>(seconds - 7 * 3600);
        const std::tm* t = std::gmtime(&local);
        price += static_cast<float>(static_cast<int>(next() % 21) - 10) / 100;
        temperature += static_cast<int>(next() % 5) - 2;
        count += 1 + static_cast<int>(next() % 3);
        const auto millis = static_cast<long long>(seconds) * 1000 + next() % 1000;
        const int quantity = 1 + static_cast<int>(next() % 50);
        // The total of the price as printed, to the cent
        const long long cents = (static_cast<long long>(std::llround(static_cast<double>(price) * 1e6)) * quantity + 5000) / 10000;
        std::snprintf(line, sizeof line, "%lld,%04d-%02d-%02d %02d:%02d:%02d,%.6f,%d.%d,%d,%lld,%d,%lld.%02lld,%d\n", static_cast<long long>(seconds),
                      t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec, static_cast<double>(price),
                      temperature / 10, temperature % 10, count, millis, quantity, cents / 100, cents % 100, (i + 3) * (i + 3));
        // Now and then a missing count of seconds
        out += i % 37 == 20 ? "NA" + std::string(std::strchr(line, ',')) : std::string(line);
    }
    return out;
}

std::string json_records() {
    std::string out = "[\n";
    char line[200];
    double lon = -122.123456, depth = 10.5;
    for (int i = 0; i < 90; ++i) {
        lon += static_cast<double>(static_cast<int>(next() % 2001) - 1000) / 1e6;
        depth = static_cast<double>(next() % 300) / 10;
        const double mag = static_cast<double>(next() % 500) / 100;
        const double lat = 37.1234567 + static_cast<double>(next() % 100000) / 1e7;
        std::snprintf(line, sizeof line, "\"mag\":%.2f,\"lon\":%.6f,\"lat\":%.7f,\"depth\":%.1f,\"title\":\"M %.2f\"},\n", mag, lon,
                      lat, depth, mag);
        // Now and then a field more at the start, which moves the place of every number after it
        out += "{\"id\":" + std::to_string(7000 + i * 3) + (i % 7 == 3 ? ",\"rank\":" + std::to_string(i % 5) : std::string()) + "," + line;
    }
    return out + "]\n";
}

std::string one_line_json() {
    std::string out = "[";
    char item[120];
    for (int i = 0; i < 150; ++i) {
        const double value = static_cast<double>(next() % 10000) / 100;
        std::snprintf(item, sizeof item, "%s{\"n\":%d,\"v\":%.2f,\"w\":[%u,%u]}", i ? "," : "", i, value, next() % 100, next() % 100);
        out += item;
    }
    return out + "]";
}

std::string log_lines() {
    std::string out;
    char line[200];
    int minute = 0, second = 0, bytes = 1000;
    for (int i = 0; i < 120; ++i) {
        second += 1 + static_cast<int>(next() % 20);
        minute += second / 60;
        second %= 60;
        bytes += static_cast<int>(next() % 4000);
        // One call per statement: the order of evaluation of arguments differs between compilers
        const uint32_t subnet = next() % 4, host = next() % 256, item = next() % 100, delta = next() % 9, id = next(), flags = next() % 4096;
        std::snprintf(line, sizeof line, "2005.06.03 22:%02d:%02d 192.168.%u.%u GET /item/%u %d -%u.5 uid=C%xq%uZ flags=0x%03x\n", minute % 60,
                      second, subnet, host, item, bytes, delta, id % 4096, id % 977, flags);
        out += line;
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
    save(dir + "numbers.csv", text(time_series()));
    save(dir + "numbers.json", text(json_records()));
    save(dir + "numbers.log", text(log_lines()));
    save(dir + "oneline.json", text(one_line_json()));

    // Every byte a marker could be, among numbers
    std::vector<uint8_t> crowded = text(log_lines());
    for (int b = 0; b < 256; ++b) {
        if (b != '\n') {
            crowded.push_back(static_cast<uint8_t>(b));
        }
    }
    save(dir + "crowded.txt", crowded);

    std::vector<uint8_t> shapes;
    double x = -122.4, y = 37.7;
    for (int s = 0; s < 120; ++s) {
        const uint32_t points = 2 + next() % 9;
        put(shapes, s, 4, false);
        put(shapes, points, 4, false);
        for (uint32_t p = 0; p < points; ++p) {
            x += static_cast<double>(static_cast<int>(next() % 2001) - 1000) * 1e-6;
            y += static_cast<double>(static_cast<int>(next() % 2001) - 1000) * 1e-6;
            put(shapes, double_bits(std::round(x * 1e6) / 1e6), 8, false);
            put(shapes, double_bits(std::round(y * 1e6) / 1e6), 8, false);
        }
    }
    save(dir + "polylines.bin", shapes);

    std::vector<uint8_t> table;
    for (int i = 0; i < 300; ++i) {
        table.push_back(static_cast<uint8_t>(next()));
    }
    int columns[3] = {};
    for (int r = 0; r < 300; ++r) {
        for (int j = 0; j < 3; ++j) {
            columns[j] += static_cast<int>(next() % 41) - 20;
            put(table, double_bits(static_cast<double>(1000 * (j + 1) + columns[j]) * 0.01), 8, false);
        }
    }
    save(dir + "header_records.bin", table);
    return 0;
}
