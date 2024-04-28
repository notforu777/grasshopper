#include <iostream>
#include <string>

const char *hex_to_bytes(const std::string &input) {
    const char *hex_chars = input.c_str();
    size_t len = input.length() / 2;
    char *output = new char[len];

    for (int i = 0; i < len; ++i) {
        std::sscanf(hex_chars + 2 * i, "%2hhx", &output[i]);
    }

    std::cout << input.length() << '\n';
    return output;
}
